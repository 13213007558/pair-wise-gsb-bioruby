# frozen_string_literal: true
#
# = bio/sequence/minimal_absent_words.rb - minimal absent words of DNA records
#
# Copyright::   Copyright (C) 2026
#               BioRuby developers
# License::     The Ruby License
#
#
module Bio
  class Sequence
    # = Bio::Sequence::MinimalAbsentWords
    #
    # A word is a _minimal absent word_ (MAW) of a collection of DNA records
    # when it never occurs in any record while every one of its proper
    # substrings does occur.  This class builds a suffix automaton over all
    # supplied records and enumerates minimal absent words in lexicographic
    # order.
    #
    # Records are never stitched together for matching: maximal ACGT runs
    # are indexed separately, so a word spanning two records, or the two
    # sides of a non-ACGT character, is never considered present.  Any
    # character outside the alphabet (for example +N+) breaks matching and
    # is never treated as a wildcard.  Occurrences may be provided by
    # different records, and each result carries verifiable occurrences of
    # both its prefix and suffix.
    #
    # Enumeration is lazy: taking the first few words with Enumerator#first
    # or breaking early from #each does not generate the whole result set.
    #
    # = Usage
    #   idx = Bio::Sequence::MinimalAbsentWords.new(['acgt', 'tgca'])
    #   idx.each(min_length: 2, max_length: 10).first(5)
    class MinimalAbsentWords
      include Enumerable

      # Canonical DNA alphabet used for matching and lexicographic ordering.
      DEFAULT_ALPHABET = 'ACGT'

      # An occurrence of a (sub)word inside one of the indexed records.
      # Coordinates are zero-based and cannot span records or separators:
      # +record.to_s[start...end_pos]+ equals the occurring word.
      Occurrence = Struct.new(:record, :record_id, :start, :end_pos,
                              keyword_init: true) do
        # Returns the occurring (sub)string.
        # ---
        # *Returns*:: String
        def word
          record.to_s[start...end_pos]
        end
      end

      # A minimal absent word with occurrences witnessing its minimality.
      #
      # +prefix_occurrence+ locates the word without its last letter;
      # +suffix_occurrence+ locates the word without its first letter.
      # Both are +nil+ for a one-letter word because the empty string has
      # no position.
      Word = Struct.new(:word, :prefix_occurrence, :suffix_occurrence,
                        keyword_init: true) do
        # Length of the absent word.
        # ---
        # *Returns*:: Integer
        def length
          word.length
        end
        alias size length
      end

      # Builds the index from one or more DNA records.
      #
      # Records may be Strings or Bio::Sequence objects.  Matching is
      # case-insensitive.  Maximal runs of alphabet characters are indexed
      # independently; every other character (including +N+, +U+ and
      # ambiguity codes) splits a run and is never a wildcard.
      # ---
      # *Arguments*:
      # * (required) _records_: one or more Strings/Bio::Sequence objects,
      #   or a single Array of them
      # * (optional keyword) _alphabet_: ordered alphabet String
      #   (default +"ACGT"+)
      # *Returns*:: Bio::Sequence::MinimalAbsentWords
      def initialize(*records, alphabet: DEFAULT_ALPHABET)
        records = records.flatten
        raise ArgumentError, 'at least one record is required' if records.empty?
        unless alphabet.is_a?(String) && alphabet.length.positive? &&
               alphabet.chars.uniq.length == alphabet.length
          raise ArgumentError,
                'alphabet must be a non-empty String of unique characters'
        end
        @records = records
        @alphabet = alphabet.upcase.freeze
        @record_ids = records.map do |rec|
          rec.entry_id if rec.respond_to?(:entry_id)
        end
        # Each segment is [record index, start in original record,
        # end position (exclusive) in original record].
        @segments = []
        records.each_with_index do |rec, rec_index|
          collect_segments(String(rec).upcase, rec_index)
        end
        build_suffix_automaton
      end

      # Records used to build the index, in the order they were given.
      # ---
      # *Returns*:: Array of String and/or Bio::Sequence objects
      attr_reader :records

      # Ordered alphabet used for matching.
      # ---
      # *Returns*:: String
      attr_reader :alphabet

      # Enumerates minimal absent words in lexicographic order, restricted
      # to _min_length_.._max_length_.  Without a block, an Enumerator is
      # returned and results are generated on demand; stopping iteration
      # early never scans longer branches or builds the whole output set.
      # ---
      # *Arguments*:
      # * (optional keyword) _min_length_: Integer (default 1)
      # * (optional keyword) _max_length_: Integer or nil (default nil)
      # *Yields*:: Bio::Sequence::MinimalAbsentWords::Word
      # *Returns*:: Enumerator when no block is given
      def each(min_length: 1, max_length: nil, &block)
        return to_enum(:each, min_length: min_length,
                               max_length: max_length) unless block_given?
        min_length = check_length(min_length, :min_length)
        max_length = check_max_length(max_length, min_length)
        enumerate(min_length, max_length, &block)
      end

      # Whether _word_ occurs in at least one indexed record.  Matching is
      # case-insensitive and cannot cross records or non-ACGT separators.
      # ---
      # *Arguments*:
      # * (required) _word_: String
      # *Returns*:: true or false
      def present?(word)
        return false if word.nil? || word.empty?
        traverse(String(word).upcase) != -1
      end

      # Whether _word_ is a minimal absent word, optionally restricted to a
      # length range.  Matching is case-insensitive.
      # ---
      # *Arguments*:
      # * (required) _word_: String
      # * (optional keyword) _min_length_: Integer (default 1)
      # * (optional keyword) _max_length_: Integer or nil (default nil)
      # *Returns*:: true or false
      def minimal_absent?(word, min_length: 1, max_length: nil)
        str = String(word).upcase
        length = str.length
        return false if length < min_length
        return false unless max_length.nil? || length <= max_length
        return false unless str.chars.all? { |ch| @alphabet.include?(ch) }
        return false if present?(str)
        prefix_ok = length == 1 || present?(str[0...-1])
        suffix_ok = length == 1 || present?(str[1..-1])
        prefix_ok && suffix_ok
      end

      private

      def check_length(value, name)
        value = Integer(value)
        raise ArgumentError, "#{name} must be positive" unless value.positive?
        value
      end

      def check_max_length(max_length, min_length)
        return nil if max_length.nil?
        max_length = check_length(max_length, :max_length)
        unless max_length >= min_length
          raise ArgumentError, 'max_length must be >= min_length'
        end
        max_length
      end

      def collect_segments(upcase_string, record_index)
        start = nil
        upcase_string.each_char.with_index do |ch, idx|
          if @alphabet.include?(ch)
            start = idx if start.nil?
          elsif start
            @segments << [record_index, start, idx]
            start = nil
          end
        end
        unless start.nil?
          @segments << [record_index, start, upcase_string.length]
        end
      end

      # Builds one text in which segments are separated by a unique byte not
      # present in any record, so transitions never create cross-record or
      # cross-N ACGT factors.
      def build_suffix_automaton
        separator = unused_separator_byte
        text = String.new.force_encoding(Encoding::ASCII_8BIT)
        # segment_start[i] is the text offset where @segments[i] begins.
        @segment_start = []
        @segments.each_with_index do |(rec_index, seg_start, seg_end), idx|
          text << separator unless idx.zero?
          @segment_start << text.length
          piece = String(@records[rec_index]).upcase[seg_start...seg_end]
          text << piece.dup.force_encoding(Encoding::ASCII_8BIT)
        end
        @length = [0]
        @link = [-1]
        @next = [{}]
        # Representative end offset of a state inside +text+.
        @end_pos = [nil]
        cur = 0
        text.each_char.with_index do |ch, idx|
          cur = sam_extend(cur, ch, idx)
        end
      end

      def unused_separator_byte
        used = Array.new(256, false)
        @records.each do |rec|
          String(rec).each_byte { |byte| used[byte] = true }
        end
        byte = used.index(false)
        raise 'no unused byte available for a record separator' if byte.nil?
        byte.chr.force_encoding(Encoding::ASCII_8BIT)
      end

      def sam_extend(cur, ch, pos)
        nxt = @length.size
        @length << @length[cur] + 1
        @link << 0
        @next << {}
        @end_pos << pos
        p = cur
        while p != -1 && !@next[p].key?(ch)
          @next[p][ch] = nxt
          p = @link[p]
        end
        if p == -1
          @link[nxt] = 0
        else
          q = @next[p][ch]
          if @length[p] + 1 == @length[q]
            @link[nxt] = q
          else
            clone = clone_state(p, q, ch)
            @link[q] = @link[nxt] = clone
          end
        end
        nxt
      end

      def clone_state(p_state, q_state, ch)
        clone = @length.size
        @length << @length[p_state] + 1
        @next << @next[q_state].dup
        @link << @link[q_state]
        @end_pos << @end_pos[q_state]
        while p_state != -1 && @next[p_state][ch] == q_state
          @next[p_state][ch] = clone
          p_state = @link[p_state]
        end
        clone
      end

      def traverse(word)
        state = 0
        word.each_char do |ch|
          state = @next[state][ch]
          return -1 if state.nil?
        end
        state
      end

      def shortest_length(state)
        state.zero? ? 0 : @length[@link[state]] + 1
      end

      def enumerate(min_length, max_length)
        # Iterative lexicographic traversal over present factors.  Frames are
        # [state, letter entering the node, next letter index] sharing one
        # mutable prefix buffer, keeping descent to depth k O(k) instead of
        # copying prefixes onto every frame.  At each letter, missing-edge
        # emissions are interleaved with child descent in alphabetical order.
        letters = @alphabet.chars
        stack = [[0, nil, 0]]
        buffer = String.new
        until stack.empty?
          state, letter, letter_index = stack.last
          buffer << letter if letter_index.zero? && !state.zero?
          depth = buffer.length
          if letter_index >= letters.length
            buffer.chop! unless state.zero?
            stack.pop
            next
          end
          current = letters[letter_index]
          stack.last[2] += 1
          child = @next[state][current]
          if child.nil?
            word = missing_edge_word(state, buffer, current, min_length,
                                     max_length)
            yield word unless word.nil?
          elsif max_length.nil? || depth < max_length
            stack << [child, current, 0]
          end
        end
      end

      def missing_edge_word(state, prefix, letter, min_length, max_length)
        word_length = prefix.length + 1
        return if word_length < min_length
        return unless max_length.nil? || word_length <= max_length
        return unless prefix.length == shortest_length(state)
        return unless state.zero? ||
                      @next[@link[state]]&.key?(letter)
        build_word(state, prefix, letter)
      end

      def build_word(state, prefix, letter)
        word_text = prefix + letter
        if state.zero?
          Word.new(word: word_text, prefix_occurrence: nil,
                           suffix_occurrence: nil)
        else
          Word.new(word: word_text,
                   prefix_occurrence: occurrence_at(state, prefix.length),
                   suffix_occurrence: occurrence_at(
                     @next[@link[state]][letter], word_text.length - 1
                   ))
        end
      end

      # Builds an Occurrence for a length-+length+ representative of +state+
      # from the representative end position stored during construction.
      def occurrence_at(state, length)
        end_offset = @end_pos[state]
        start_offset = end_offset - length + 1
        segment_index = segment_index_at(start_offset)
        rec_index, seg_start, = @segments[segment_index]
        base = @segment_start[segment_index]
        Occurrence.new(record: @records[rec_index],
                       record_id: @record_ids[rec_index],
                       start: seg_start + (start_offset - base),
                       end_pos: seg_start + (end_offset - base) + 1)
      end

      # Segment containing a text offset.  Offsets of emitted witnesses
      # always lie inside an ACGT segment.
      def segment_index_at(offset)
        index = (@segment_start.bsearch_index { |start| start > offset } ||
                 @segment_start.length) - 1
        raise 'offset outside indexed segments' if index.negative?
        index
      end
    end
  end
end
