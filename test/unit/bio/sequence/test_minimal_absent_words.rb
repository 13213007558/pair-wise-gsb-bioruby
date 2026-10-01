# frozen_string_literal: true

# loading helper routine for testing bioruby
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'test/unit'
require 'bio/sequence'
require 'bio/sequence/na'
require 'bio/sequence/minimal_absent_words'

module Bio
  class TestSequenceMinimalAbsentWordsIndex < Test::Unit::TestCase
    def test_enumerates_minimal_absent_words_in_order
      index = Bio::Sequence::MinimalAbsentWordsIndex.new('atgc')

      assert_equal(%w[aa ac ag ca cc cg ct ga gg gt ta tc tt],
                   index.map(&:word))
    end

    def test_filters_by_length_range
      index = Bio::Sequence::MinimalAbsentWordsIndex.new('atgc')

      assert_equal([], index.each(3, 4).to_a)
      assert_equal(%w[aa ac ag ca cc cg ct ga gg gt ta tc tt],
                   index.each(1, 2).map(&:word))
    end

    def test_brute_force_comparison_for_small_records
      records = ['acNgt', 'ca']
      index = Bio::Sequence::MinimalAbsentWordsIndex.new(records)

      assert_equal(brute_force_minimal_absent_words(records, 1, 4),
                   index.each(1, 4).map(&:word))
      index.each(1, 4) do |word|
        assert_evidence(records, word)
      end
    end

    def test_non_acgt_characters_split_matching_runs
      combined = Bio::Sequence::MinimalAbsentWordsIndex.new('acNgt')
      separated = Bio::Sequence::MinimalAbsentWordsIndex.new('ac', 'gt')

      assert_equal(separated.map(&:word), combined.map(&:word))
      assert_include(combined.map(&:word), 'cg')
      assert_include(combined.map(&:word), 'at')
    end

    def test_records_are_not_joined
      index = Bio::Sequence::MinimalAbsentWordsIndex.new('ac', 'gt')

      assert_include(index.map(&:word), 'cg')
      assert_not_include(index.each(1, 4).map(&:word), 'acgt')
    end

    def test_evidence_can_be_checked_against_input_records
      word = Bio::Sequence::MinimalAbsentWordsIndex
             .new('acNgt', 'ca')
             .find { |candidate| candidate.word == 'aca' }

      assert_equal('ac', occurrence_text('acNgt', word.prefix_occurrence, 2))
      assert_equal('ca', occurrence_text('ca', word.suffix_occurrence, 2))
    end

    def test_supports_prefix_and_suffix_from_different_records
      records = ['gt', 'ta']
      word = Bio::Sequence::MinimalAbsentWordsIndex
             .new(records)
             .find { |candidate| candidate.word == 'gta' }

      assert_equal(0, word.prefix_occurrence.record_index)
      assert_equal(1, word.suffix_occurrence.record_index)
      assert_equal('gt', occurrence_text(records[0],
                                          word.prefix_occurrence, 2))
      assert_equal('ta', occurrence_text(records[1],
                                          word.suffix_occurrence, 2))
    end

    def test_single_letter_words_have_no_proper_substring_evidence
      word = Bio::Sequence::MinimalAbsentWordsIndex.new('atg').find do |candidate|
        candidate.word == 'c'
      end

      assert_equal('c', word.word)
      assert_nil(word.prefix_occurrence)
      assert_nil(word.suffix_occurrence)
    end

    def test_returns_enumerator_without_block_and_is_lazy
      enumerator = Bio::Sequence::MinimalAbsentWordsIndex
                   .new('a' + 'c' * 1000)
                   .each(1, nil)

      assert_instance_of(Enumerator, enumerator)
      assert_equal('aa', enumerator.first.word)
    end

    def test_long_sequence_reaches_large_length_without_walking_all_candidates
      sequence = ('a' * 5000) + 'cg'
      enumerator = Bio::Sequence::MinimalAbsentWordsIndex
                   .new(sequence)
                   .each(5001, 5001)

      assert_equal('a' * 5001, enumerator.first.word)
    end

    def test_na_convenience_api
      sequence = Bio::Sequence::NA.new('ATGC')

      assert_equal(%w[aa ac ag ca cc cg ct ga gg gt ta tc tt],
                   sequence.minimal_absent_words.map(&:word))
    end

    def test_biosequence_convenience_api
      sequence = Bio::Sequence.new('ATGC')

      assert_equal(%w[aa ac ag ca cc cg ct ga gg gt ta tc tt],
                   sequence.minimal_absent_words.map(&:word))
    end

    def test_validates_length_arguments
      index = Bio::Sequence::MinimalAbsentWordsIndex.new('ac')

      assert_raise(ArgumentError) { index.each(0).to_a }
      assert_raise(ArgumentError) { index.each(3, 2).to_a }
    end

    private

    def occurrence_text(record, occurrence, length)
      record[occurrence.offset, length]
    end

    def assert_evidence(records, word)
      return if word.word.length == 1

      prefix = word.word[0...-1]
      suffix = word.word[1..-1]
      prefix_occurrence = word.prefix_occurrence
      suffix_occurrence = word.suffix_occurrence

      assert_equal(prefix,
                   occurrence_text(records[prefix_occurrence.record_index],
                                   prefix_occurrence, prefix.length))
      assert_equal(suffix,
                   occurrence_text(records[suffix_occurrence.record_index],
                                   suffix_occurrence, suffix.length))
    end

    def brute_force_minimal_absent_words(records, min_length, max_length)
      factors = Set.new
      records.each do |record|
        record.downcase.delete(" \t\n\r").split(/[^acgt]/).each do |segment|
          (0...segment.length).each do |start|
            (start + 1..segment.length).each do |ending|
              factors << segment[start...ending]
            end
          end
        end
      end

      (min_length..max_length).flat_map do |length|
        %w[a c g t].repeated_permutation(length).map(&:join).select do |word|
          !factors.include?(word) &&
            (length == 1 ||
             (factors.include?(word[0...-1]) &&
              factors.include?(word[1..-1])))
        end
      end.sort
    end
  end
end
