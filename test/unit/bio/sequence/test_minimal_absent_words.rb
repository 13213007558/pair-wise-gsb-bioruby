#
# = test/unit/bio/sequence/test_minimal_absent_words.rb
#
# Copyright::   Copyright (C) 2026
#               BioRuby developers
# License::     The Ruby License
#

# loading helper routine for testing bioruby
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'test/unit'
require 'bio/sequence'
require 'bio/sequence/na'
require 'bio/sequence/minimal_absent_words'

module Bio
  class TestSequenceMinimalAbsentWords < Test::Unit::TestCase
    ALPHABET = %w[A C G T].freeze

    def maw_bruteforce(records, min_length, max_length)
      factors = { '' => true }
      records.each do |record|
        record.to_s.upcase.split(/[^ACGT]/).each do |segment|
          (0...segment.length).each do |i|
            (i...segment.length).each do |j|
              factors[segment[i..j]] = true
            end
          end
        end
      end
      result = []
      (min_length..max_length).each do |length|
        words = if length == 1
                  ALPHABET
                else
                  ALPHABET.product(*Array.new(length - 1) { ALPHABET })
                end
        words.each do |letters|
          word = letters.is_a?(Array) ? letters.join : letters
          next if factors[word]
          prefix_ok = length == 1 || factors[word[0...-1]]
          suffix_ok = length == 1 || factors[word[1..-1]]
          result << word if prefix_ok && suffix_ok
        end
      end
      result.sort
    end

    def words_of(index, min_length: 1, max_length: 6)
      index.each(min_length: min_length, max_length: max_length).map(&:word)
    end

    def test_acgt
      index = Bio::Sequence::MinimalAbsentWords.new('acgt')
      expected = %w[AA AG AT CA CC CT GA GC GG TA TC TG TT]
      assert_equal(expected, words_of(index))
    end

    def test_lexicographic_order
      index = Bio::Sequence::MinimalAbsentWords.new('ttggccaa')
      words = index.each(max_length: 4).map(&:word)
      assert_equal(words.sort, words)
    end

    def test_length_range
      index = Bio::Sequence::MinimalAbsentWords.new('acgt')
      assert_equal([], index.each(min_length: 1, max_length: 1).to_a)
      words = index.each(min_length: 2, max_length: 2).map(&:word)
      assert_equal(13, words.length)
      assert(words.all? { |word| word.length == 2 })
      assert_equal(%w[AA AG AT], index.each.first(3).map(&:word))
    end

    def test_records_are_not_concatenated
      # 'ac' and 'gt' together must not create 'cg'.
      index = Bio::Sequence::MinimalAbsentWords.new(['ac', 'gt'])
      assert_equal(false, index.present?('cg'))
      assert(index.each.map(&:word).include?('CG'))
      index_single = Bio::Sequence::MinimalAbsentWords.new('acgt')
      assert_equal(true, index_single.present?('cg'))
      assert_equal(false, index_single.each.map(&:word).include?('CG'))
    end

    def test_non_acgt_characters_break_matching
      # N is a separator, not a wildcard: neither 'an' style matches nor
      # words crossing the N are present.
      index = Bio::Sequence::MinimalAbsentWords.new('acgNtgca')
      assert_equal(true, index.present?('tgca'))
      assert_equal(true, index.present?('acg'))
      assert_equal(false, index.present?('acgnt'))
      assert_equal(false, index.present?('gt'))
      # As separated segments 'acg' and 'tgca', 'GA' is minimal absent.
      assert(index.each(max_length: 2).map(&:word).include?('GA'))
    end

    def test_ambiguous_characters_are_separators
      index = Bio::Sequence::MinimalAbsentWords.new('aRmkswy')
      # Only 'A' is present; it is a single-letter segment, so AA is also
      # absent.  Every ambiguity code cuts matching.
      assert_equal(%w[AA C G T], index.each.map(&:word))
    end

    def test_no_acgt_content
      index = Bio::Sequence::MinimalAbsentWords.new(['NNN', ''])
      assert_equal(%w[A C G T], index.each.map(&:word))
      assert_equal(false, index.present?('A'))
    end

    def test_witness_occurrences_are_verifiable
      index = Bio::Sequence::MinimalAbsentWords.new(['ac', 'gt'])
      word = index.each.find { |item| item.word == 'CT' }
      records = index.records
      assert_equal('c', word.prefix_occurrence.word)
      assert_equal('t', word.suffix_occurrence.word)
      # Prefix and suffix can come from different records.
      assert_same(records[0], word.prefix_occurrence.record)
      assert_same(records[1], word.suffix_occurrence.record)
      assert_equal([1, 2], [word.prefix_occurrence.start,
                           word.prefix_occurrence.end_pos])
      assert_equal([1, 2], [word.suffix_occurrence.start,
                           word.suffix_occurrence.end_pos])
      # Zero-based coordinates slice exactly the witnessed factors.
      prefix = word.prefix_occurrence
      suffix = word.suffix_occurrence
      assert_equal('c', prefix.record[prefix.start...prefix.end_pos])
      assert_equal('t', suffix.record[suffix.start...suffix.end_pos])
    end

    def test_witness_coordinates_inside_record_with_separator
      index = Bio::Sequence::MinimalAbsentWords.new('acNtg')
      index.each(max_length: 3) do |word|
        next if word.length == 1
        [word.prefix_occurrence, word.suffix_occurrence].each do |occ|
          assert_equal(occ.word.upcase,
                       'acNtg'[occ.start...occ.end_pos].upcase)
          assert(!occ.word.upcase.include?('N'))
        end
      end
    end

    def test_one_letter_words_have_no_occurrences
      word = Bio::Sequence::MinimalAbsentWords.new('').each.first
      assert_equal('A', word.word)
      assert_nil(word.prefix_occurrence)
      assert_nil(word.suffix_occurrence)
    end

    def test_record_id_is_exposed
      sequence = Bio::Sequence.new('acgt')
      sequence.entry_id = 'rec_1'
      index = Bio::Sequence::MinimalAbsentWords.new(sequence)
      word = index.each.find { |item| item.word == 'AA' }
      assert_equal('rec_1', word.prefix_occurrence.record_id)
    end

    def test_present_and_minimal_absent_predicates
      index = Bio::Sequence::MinimalAbsentWords.new('acgt')
      assert_equal(true, index.present?('Ac'))
      assert_equal(false, index.present?('aa'))
      assert_equal(true, index.minimal_absent?('aa'))
      assert_equal(false, index.minimal_absent?('aaa'))
      assert_equal(false, index.minimal_absent?('ac'))
      assert_equal(false, index.minimal_absent?('an'))
      assert_equal(true, index.minimal_absent?('aa', min_length: 2))
      assert_equal(false, index.minimal_absent?('aa', min_length: 3))
    end

    def test_enumerator_is_lazy
      sequence = 'A' * 10_000
      index = Bio::Sequence::MinimalAbsentWords.new(sequence)
      yielded = 0
      index.each do |word|
        yielded += 1
        assert_equal(10_001, word.length)
        break
      end
      assert_equal(1, yielded)
      assert_equal(10_001, index.each.first.length)
      # The first MAW must be found without materializing 4^k candidates;
      # the whole index construction + first lookup stays well under a
      # brute-force enumeration budget.
      assert_operator(index.each.first.word.length, :==, 10_001)
    end

    def test_block_enumeration_matches_enumerator
      index = Bio::Sequence::MinimalAbsentWords.new('acgtacgt')
      collected = []
      index.each(max_length: 4) { |word| collected << word.word }
      assert_equal(index.each(max_length: 4).map(&:word), collected)
    end

    def test_na_shortcut
      na = Bio::Sequence::NA.new('acgt')
      assert_equal(Bio::Sequence::MinimalAbsentWords.new(na).each.map(&:word),
                   na.minimal_absent_words.map(&:word))
      assert_equal(Bio::Sequence::MinimalAbsentWords.new(na).each.map(&:word),
                   na.maw.map(&:word))
    end

    def test_argument_validation
      assert_raise(ArgumentError) { Bio::Sequence::MinimalAbsentWords.new }
      assert_raise(ArgumentError) do
        Bio::Sequence::MinimalAbsentWords.new('a', alphabet: 'AA')
      end
      assert_raise(ArgumentError) do
        Bio::Sequence::MinimalAbsentWords.new('a').each(min_length: 0).to_a
      end
      assert_raise(ArgumentError) do
        Bio::Sequence::MinimalAbsentWords.new('a')
          .each(min_length: 3, max_length: 2).to_a
      end
    end

    def test_random_collections_match_bruteforce
      srand(1234)
      100.times do
        records = Array.new(rand(1..3)) do
          length = rand(0..14)
          Array.new(length) do
            rand(10) < 7 ? ALPHABET.sample : %w[N U -].sample
          end.join
        end
        index = Bio::Sequence::MinimalAbsentWords.new(records)
        words = index.each(min_length: 1, max_length: 5).map(&:word)
        assert_equal(maw_bruteforce(records, 1, 5), words)
        assert_equal(words.sort, words)
      end
    end

    def test_random_witnesses_match
      srand(4321)
      30.times do
        records = Array.new(2) do
          length = rand(2..12)
          Array.new(length) { ALPHABET.sample }.join
        end
        index = Bio::Sequence::MinimalAbsentWords.new(records)
        index.each(max_length: 4) do |word|
          next if word.word.length == 1
          prefix = word.prefix_occurrence
          suffix = word.suffix_occurrence
          assert_equal(word.word[0...-1], prefix.word.upcase)
          assert_equal(word.word[1..-1], suffix.word.upcase)
          assert(prefix.record.to_s[prefix.start...prefix.end_pos]
                   .upcase.index(/[^ACGT]/).nil?)
          assert(suffix.record.to_s[suffix.start...suffix.end_pos]
                   .upcase.index(/[^ACGT]/).nil?)
        end
      end
    end
  end
end
