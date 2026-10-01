# frozen_string_literal: true

require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 2,
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'test/unit'
require 'bio/codon_optimizer'
require 'bio/sequence/na'

module Bio
  module CodonOptimizerTestHelper
    A_WEIGHTS = {
      'gct' => 10,
      'gcc' => 9,
      'gca' => 8,
      'gcg' => 7
    }.freeze

    F_WEIGHTS = {
      'ttt' => 10,
      'ttc' => 1
    }.freeze

    module_function

    def optimize(sequence, weights, forbidden, **options)
      Bio::CodonOptimizer.optimize(
        sequence, 1,
        **options,
        codon_weights: weights,
        forbidden_substrings: forbidden
      )
    end

    def best_bruteforce(candidates, weights, forbidden, gc_range)
      feasible = candidates.select do |sequence|
        forbidden.none? { |pattern| sequence.include?(pattern) } &&
          gc_range.cover?(sequence.count('gc'))
      end
      best = feasible.max_by { |sequence| score(sequence, weights) }
      { sequence: best, score: score(best, weights) }
    end

    def score(sequence, weights)
      sequence.scan(/.{3}/).sum { |codon| weights.fetch(codon, 0) }
    end
  end

  class TestCodonOptimizer < Test::Unit::TestCase
    include CodonOptimizerTestHelper

    def test_exhaustive_short_peptide_optimum
      candidates = %w[gct gcc gca gcg].product(%w[ttt ttc]).map(&:join)
      forbidden = %w[tttt ctt]
      weights = A_WEIGHTS.merge(F_WEIGHTS)
      optimum = best_bruteforce(candidates, weights, forbidden, 2..5)

      result = Bio::CodonOptimizer.optimize(
        'gctttt', 1,
        codon_weights: weights,
        forbidden_substrings: forbidden,
        gc_count: 2..5
      )

      assert(result.feasible?)
      assert_equal(optimum.fetch(:sequence), result.sequence)
      assert_equal(optimum.fetch(:score), result.score)
      assert_equal('AF', result.sequence.translate(1, 1))
    end

    def test_cross_codon_forbidden_substring_changes_earlier_choice
      weights = A_WEIGHTS.merge(F_WEIGHTS)
      result = optimize('gctttt', weights, ['cttt'])

      assert(result.feasible?)
      assert_equal('gcattt', result.sequence)
      assert_equal(18, result.score)
      assert_equal('AF', result.sequence.translate(1, 1))
    end

    def test_fixed_position_cannot_change_even_with_better_weight
      result = optimize('gctttt', A_WEIGHTS.merge(F_WEIGHTS), [],
                        fixed_positions: [0])

      assert(result.feasible?)
      assert_equal('gctttt', result.sequence)
      assert_equal(20, result.score)
    end

    def test_range_fixed_positions_and_string_pattern
      result = Bio::CodonOptimizer.optimize(
        'gctttt', 1,
        codon_weights: A_WEIGHTS.merge(F_WEIGHTS),
        forbidden_substrings: 'cttt',
        fixed_positions: 0..1
      )

      assert_equal(false, result.feasible?)
      assert_equal('forbidden_substring', result.conflicts.dig(0, 'type'))
    end

    def test_fixed_position_can_make_forbidden_constraint_infeasible
      result = optimize('gctttt', A_WEIGHTS.merge(F_WEIGHTS), ['cttt'],
                        fixed_positions: [0, 1])

      assert_equal(false, result.feasible?)
      assert_nil(result.sequence)
      assert_equal('forbidden_substring', result.conflicts.dig(0, 'type'))
    end

    def test_gc_conflict_reports_achievable_range
      result = optimize('tttttt', F_WEIGHTS, [], gc_count: 3..4)

      assert_equal(false, result.feasible?)
      assert_equal('gc_count', result.conflicts.dig(0, 'type'))
      assert_equal(3, result.conflicts.dig(0, 'requested_min'))
      assert_equal(4, result.conflicts.dig(0, 'requested_max'))
      assert_equal(0, result.conflicts.dig(0, 'achievable_min'))
      assert_equal(2, result.conflicts.dig(0, 'achievable_max'))
    end

    def test_gc_range_alias
      result = optimize('gctttt', A_WEIGHTS.merge(F_WEIGHTS), [],
                        gc_range: 2..3)

      assert(result.feasible?)
      assert_equal(2, result.gc_count)
    end

    def test_caller_supplied_genetic_code_table
      table = Bio::CodonTable[2]
      result = Bio::CodonOptimizer.optimize(
        'atgtga', table,
        codon_weights: { 'atg' => 100, 'tga' => 1, 'tgg' => 10 },
        start_codon_strategy: :optimize
      )

      assert(result.feasible?)
      assert_equal('atgtgg', result.sequence)
      assert_equal('MW', result.sequence.translate(1, table))
    end

    def test_start_and_stop_codons_are_preserved_by_default
      result = Bio::CodonOptimizer.optimize(
        'gtggcttag', 11,
        codon_weights: { 'gtt' => 100, 'taa' => 100, 'gct' => 1 },
        forbidden_substrings: []
      )

      assert(result.feasible?)
      assert_equal('gtggcttag', result.sequence)
      assert_equal('VA*', result.sequence.translate(1, 11))
    end

    def test_start_and_stop_codons_can_be_optimized_when_requested
      result = Bio::CodonOptimizer.optimize(
        'gtggcttag', 11,
        codon_weights: { 'gtt' => 100, 'taa' => 100, 'gct' => 1 },
        start_codon_strategy: :optimize,
        stop_codon_strategy: :optimize
      )

      assert(result.feasible?)
      assert_equal('gttgcttaa', result.sequence)
      assert_equal('VA*', result.sequence.translate(1, 11))
    end

    def test_rna_input_returns_rna_and_retranslates_with_caller_table
      result = Bio::CodonOptimizer.optimize(
        'GCUUUU', 1,
        codon_weights: { 'gcu' => 10, 'uuu' => 10, 'uuc' => 1 },
        forbidden_substrings: ['cuuu']
      )

      assert(result.feasible?)
      assert_equal('gcauuu', result.sequence)
      assert_equal('AF', result.sequence.translate(1, 1))
    end

    def test_sequence_convenience_method
      sequence = Bio::Sequence::NA.new('gctttt')
      result = sequence.optimize_codons(
        1,
        codon_weights: A_WEIGHTS.merge(F_WEIGHTS),
        forbidden_substrings: ['cttt']
      )

      assert_equal('gcattt', result.sequence)
    end

    def test_invalid_sequence_length_raises
      assert_raise(ArgumentError) do
        optimize('gcttt', F_WEIGHTS, [])
      end
    end
  end
end
