# frozen_string_literal: true

#
# = bio/sequence/minimal_absent_words.rb - minimal absent words of DNA sequences
#
# Copyright:: Copyright (C) 2026 BioRuby project
# License:: The Ruby License
#

module Bio
  require 'bio/sequence' unless const_defined?(:Sequence)

class Sequence

# Indexes one or more DNA records and enumerates their minimal absent words
# in lexicographic order. A minimal absent word is absent even though every
# proper substring of the word occurs in at least one indexed record.
#
# The DNA alphabet is limited to A, C, G, and T. Lowercase and uppercase
# letters are equivalent. Characters outside the alphabet, including N, do
# not match any nucleotide and split a record into separate matching runs.
# Records are never joined, so a substring cannot be supported across a
# record boundary.
class MinimalAbsentWordsIndex
  include Enumerable

  # DNA letters in lexicographic order after normalization.
  ALPHABET = %w[a c g t].freeze

  # A verifiable occurrence of a prefix or suffix used as evidence.
  # +record_index+ identifies the input record and +offset+ is the zero-based
  # start position in that normalized record.
  Occurrence = Struct.new(:record_index, :offset, keyword_init: true)

  # A minimal absent word and evidence for its two longest proper substrings.
  Word = Struct.new(:word,
                    :prefix_occurrence,
                    :suffix_occurrence,
                    keyword_init: true)

  # Builds the index. Each argument may be one sequence or an array of
  # sequences. Bio::Sequence, Bio::Sequence::NA, and plain String values are
  # accepted.
  def initialize(*sequences)
    @input_records = sequences.flatten.map { |sequence| normalize(sequence) }
    @transitions = [{}]
    @suffix_links = [-1]
    @min_lengths = [0]
    @max_lengths = [0]
    @first_occurrences = [nil]
    build_index
    calculate_minimum_lengths
    build_shortest_factor_trie
  end

  # Returns an Enumerator when no block is supplied. The first requested
  # result is generated lazily; the complete result set is never generated
  # merely to satisfy Enumerable methods such as +take+ or +first+.
  #
  # ---
  # *Arguments*:
  # * (optional) _min_length_: minimum word length (default 1)
  # * (optional) _max_length_: maximum word length, or nil for no maximum
  # *Yields*:: Bio::Sequence::MinimalAbsentWordsIndex::Word
  # *Returns*:: Enumerator without a block, otherwise nil
  def each(min_length = 1, max_length = nil)
    return enum_for(:each, min_length, max_length) unless block_given?

    validate_length_range(min_length, max_length)
    enumerate(min_length, max_length) { |word| yield word }
    nil
  end

  # Creates a temporary index and enumerates minimal absent words in the
  # supplied records.
  def self.each(sequences, min_length = 1, max_length = nil, &block)
    new(sequences).each(min_length, max_length, &block)
  end

  private

  def normalize(sequence)
    String(sequence).downcase.delete(" \t\n\r")
  end

  def build_index
    trie_children = [{}]
    trie_depths = [0]
    trie_occurrences = [nil]
    trie_to_sam_state = { 0 => 0 }

    @input_records.each_with_index do |record, record_index|
      trie_node = 0
      record.each_char.with_index do |character, offset|
        unless ALPHABET.include?(character)
          trie_node = 0
          next
        end

        child = trie_children[trie_node][character]
        unless child
          child = trie_children.size
          trie_children << {}
          trie_depths << trie_depths[trie_node] + 1
          trie_occurrences << [record_index, offset]
          trie_children[trie_node][character] = child
        end
        trie_node = child
      end
    end

    queue = [0]
    until queue.empty?
      parent = queue.shift
      trie_children[parent].each do |character, trie_node|
        extend_trie(character, parent, trie_node, trie_depths[trie_node],
                    trie_occurrences[trie_node], trie_to_sam_state)
        queue << trie_node
      end
    end
  end

  def extend_trie(character, trie_parent, trie_node_id, trie_depth,
                  occurrence, trie_to_sam_state)
    suffix_parent = trie_to_sam_state[trie_parent]
    existing = @transitions[suffix_parent][character]
    if existing
      if @max_lengths[existing] == trie_depth
        trie_to_sam_state[trie_node_id] = existing
      else
        clone = clone_state(existing, trie_depth)
        redirect_transitions(suffix_parent, character, existing, clone)
        @suffix_links[existing] = clone
        trie_to_sam_state[trie_node_id] = clone
      end
      return
    end

    current = add_state(trie_depth)
    @first_occurrences[current] = occurrence
    trie_to_sam_state[trie_node_id] = current

    state = suffix_parent
    while state != -1 && @transitions[state][character].nil?
      @transitions[state][character] = current
      state = @suffix_links[state]
    end

    if state == -1
      @suffix_links[current] = 0
      return
    end

    existing = @transitions[state][character]
    if @max_lengths[state] + 1 == @max_lengths[existing]
      @suffix_links[current] = existing
      return
    end

    clone = clone_state(existing, @max_lengths[state] + 1)
    redirect_transitions(state, character, existing, clone)
    @suffix_links[existing] = clone
    @suffix_links[current] = clone
  end

  def redirect_transitions(state, character, existing, clone)
    while state != -1 && @transitions[state][character] == existing
      @transitions[state][character] = clone
      state = @suffix_links[state]
    end
  end

  def add_state(max_length)
    @transitions << {}
    @suffix_links << nil
    @min_lengths << nil
    @max_lengths << max_length
    @first_occurrences << nil
    @transitions.size - 1
  end

  def clone_state(source, max_length)
    state = add_state(max_length)
    @transitions[state] = @transitions[source].dup
    @suffix_links[state] = @suffix_links[source]
    @first_occurrences[state] = @first_occurrences[source]&.dup
    state
  end

  def calculate_minimum_lengths
    (1...@transitions.size).each do |state|
      @min_lengths[state] = @max_lengths[@suffix_links[state]] + 1
    end
  end

  def build_shortest_factor_trie
    @trie_parents = Array.new(@transitions.size)
    incoming = Array.new(@transitions.size) { [] }

    @transitions.each_with_index do |transitions, source|
      transitions.each do |character, target|
        next unless character.is_a?(String)

        incoming[target] << [source, character]
      end
    end

    (1...@transitions.size).each do |state|
      prefix_length = @min_lengths[state] - 1
      parent, _character = incoming[state].find do |source, _|
        @min_lengths[source] <= prefix_length &&
          prefix_length <= @max_lengths[source]
      end
      @trie_parents[state] = parent
    end
  end

  def enumerate(min_length, max_length)
    stack = [[0, 0, 0]]

    until stack.empty?
      state, depth, letter_index = stack.last
      if letter_index == ALPHABET.size
        stack.pop
        next
      end

      stack.last[2] += 1
      character = ALPHABET[letter_index]
      target = @transitions[state][character]

      if target.nil?
        candidate_length = depth + 1
        next if candidate_length < min_length
        next if max_length && candidate_length > max_length
        next unless missing_transition_is_minimal?(state, character)

        yield build_word(state, character)
      elsif @trie_parents[target] == state &&
            @min_lengths[target] == depth + 1 &&
            (max_length.nil? || depth + 1 <= max_length)
        stack << [target, depth + 1, 0]
      end
    end
  end

  def missing_transition_is_minimal?(state, character)
    state.zero? ||
      @transitions[@suffix_links[state]].key?(character)
  end

  def build_word(state, character)
    return Word.new(word: character,
                    prefix_occurrence: nil,
                    suffix_occurrence: nil) if state.zero?

    prefix_length = @min_lengths[state]
    prefix = occurrence(state, prefix_length)
    suffix_state = @transitions[@suffix_links[state]][character]
    suffix = occurrence(suffix_state, prefix_length)

    Word.new(word: @input_records[prefix.record_index]
                       .slice(prefix.offset, prefix_length) + character,
             prefix_occurrence: prefix,
             suffix_occurrence: suffix)
  end

  def occurrence(state, length)
    record_index, ending = @first_occurrences[state]
    Occurrence.new(record_index: record_index,
                   offset: ending - length + 1)
  end

  def validate_length_range(min_length, max_length)
    unless min_length.is_a?(Integer) && min_length.positive?
      raise ArgumentError, 'min_length must be a positive Integer'
    end
    return if max_length.nil?

    unless max_length.is_a?(Integer) && max_length >= min_length
      raise ArgumentError,
            'max_length must be an Integer no smaller than min_length'
    end
  end
end

end # Sequence
end # Bio
