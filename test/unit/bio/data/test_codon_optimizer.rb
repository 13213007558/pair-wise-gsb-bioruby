#
# test/unit/bio/data/test_codon_optimizer.rb - Unit test for CodonOptimizer
#
# Copyright:: Copyright (C) 2026 BioRuby project
# License:: The Ruby License
#
# frozen_string_literal: true


require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'test/unit'
require 'bio/data/codon_optimizer'
require 'bio/sequence'
require 'bio/sequence/na'

module Bio
  class TestCodonOptimizer < Test::Unit::TestCase
    WEIGHTS = {
      'ttt' => 10, 'ttc' => 1,
      'gct' => 10, 'gcc' => 2, 'gca' => 8, 'gcg' => 9,
      'taa' => 10, 'tag' => 1, 'tga' => 5
    }.freeze

    FORBIDDEN = %w[ttg tgc ttgc tttg].freeze

    def setup
      @input = 'TTTGCTTAA'
      @choices = choices_for(@input)
      @optimizer = Bio::CodonOptimizer.new(
        codon_table: 1,
        codon_weights: WEIGHTS,
        forbidden_substrings: FORBIDDEN,
        fixed_positions: [1],
        fix_start: false,
        fix_stop: false
      )
    end

    def test_optimize_finds_enumerated_optimum
      expected = brute_force(@input, WEIGHTS, FORBIDDEN, [1])
      result = @optimizer.optimize(@input)

      assert_equal(true, result.feasible?)
      assert_equal(expected[:sequence], result.sequence)
      assert_equal(expected[:score], result.score)
      assert_equal(gc_count(expected[:sequence]), result.gc_count)
      assert_equal('FA*', translate(result.sequence))
    end

    def test_checks_forbidden_substring_across_codon_boundary
      result = Bio::CodonOptimizer.optimize(
        'TTTGCTTAA',
        codon_table: 1,
        codon_weights: WEIGHTS,
        forbidden_substrings: %w[tttgct],
        fix_start: false,
        fix_stop: false
      )

      expected = brute_force('TTTGCTTAA', WEIGHTS, %w[tttgct], [])
      assert_equal(expected[:sequence], result.sequence)
      assert_equal(expected[:score], result.score)
      assert_not_equal('tttgcttaa', result.sequence)
    end

    def test_fixed_positions_are_preserved
      result = Bio::CodonOptimizer.optimize(
        'TTTGCTTAA',
        codon_table: 1,
        codon_weights: WEIGHTS,
        fixed_positions: [1],
        fix_start: false,
        fix_stop: false
      )

      assert_equal('gct', result.sequence[3, 3])
      assert_equal('FA*', translate(result.sequence))
    end

    def test_gc_range_selects_optimal_bounded_candidate
      expected = brute_force(@input, WEIGHTS, [], [], 5..5)
      result = Bio::CodonOptimizer.optimize(
        @input,
        codon_table: 1,
        codon_weights: WEIGHTS,
        gc_range: 5..5,
        fix_start: false,
        fix_stop: false
      )

      assert_equal(expected[:sequence], result.sequence)
      assert_equal(5, result.gc_count)
    end

    def test_returns_forbidden_conflict_without_raising
      result = Bio::CodonOptimizer.optimize(
        'TTTGGTTAA',
        codon_table: 1,
        forbidden_substrings: %w[gg],
        fix_start: false,
        fix_stop: false
      )

      assert_equal(false, result.feasible?)
      assert_nil(result.sequence)
      assert_equal('forbidden_substrings_impossible',
                   result.conflicts.first[:type])
    end

    def test_returns_gc_conflict_after_constrained_search
      result = Bio::CodonOptimizer.optimize(
        'TTTGCTTAA',
        codon_table: 1,
        forbidden_substrings: FORBIDDEN,
        gc_range: 0..0,
        fix_start: false,
        fix_stop: false
      )

      assert_equal(false, result.feasible?)
      assert_equal('gc_range_impossible', result.conflicts.first[:type])
    end

    def test_rna_input_returns_rna_and_translates_after_dna_conversion
      result = @optimizer.optimize('UUUGCUUAA')

      assert_equal(true, result.feasible?)
      assert_match(/\A[acgu]+\z/, result.sequence)
      assert_equal('FA*', translate(result.sequence.tr('u', 't')))
    end

    def test_start_and_stop_are_fixed_by_default
      result = Bio::CodonOptimizer.optimize(
        'ATGGCTTAA',
        codon_table: 1,
        codon_weights: WEIGHTS
      )

      assert_equal('atg', result.sequence[0, 3])
      assert_equal('taa', result.sequence[6, 3])
    end

    def test_start_and_stop_options_are_independent
      changed_start = Bio::CodonOptimizer.optimize(
        'TTTGCTTAA',
        codon_table: 1,
        codon_weights: { 'ttc' => 10, 'tga' => 100 },
        fix_start: false,
        fix_stop: true
      )
      assert_equal('ttc', changed_start.sequence[0, 3])
      assert_equal('taa', changed_start.sequence[6, 3])

      changed_stop = Bio::CodonOptimizer.optimize(
        'ATGGCTTAA',
        codon_table: 1,
        codon_weights: { 'tga' => 100 },
        fix_start: true,
        fix_stop: false
      )
      assert_equal('atg', changed_stop.sequence[0, 3])
      assert_equal('tga', changed_stop.sequence[6, 3])
    end

    def test_stop_uses_supplied_genetic_code
      result = Bio::CodonOptimizer.optimize(
        'ATGGCTTGA',
        codon_table: 6,
        codon_weights: { 'tga' => 10 },
        fix_start: true,
        fix_stop: false
      )

      assert_equal('tga', result.sequence[6, 3])
      assert_equal('MA*', translate(result.sequence, Bio::CodonTable[6]))
    end

    def test_invalid_length_is_reported_as_conflict
      result = Bio::CodonOptimizer.optimize('atg', codon_table: 1)
      assert_equal(true, result.feasible?)

      result = Bio::CodonOptimizer.optimize('at', codon_table: 1)
      assert_equal(false, result.feasible?)
      assert_equal('invalid_sequence_length',
                   result.conflicts.first[:type])
    end

    def test_out_of_range_fixed_position_is_reported_as_conflict
      result = Bio::CodonOptimizer.optimize(
        'ATGGCTTAA',
        codon_table: 1,
        fixed_positions: [3]
      )

      assert_equal(false, result.feasible?)
      assert_equal('fixed_position_out_of_range',
                   result.conflicts.first[:type])
    end

    def test_caller_supplied_codon_table_is_accepted
      table = Bio::CodonTable.copy(1)
      result = Bio::CodonOptimizer.optimize(
        'ATGGCTTAA',
        codon_table: table,
        fix_start: true,
        fix_stop: true
      )

      assert_equal('atggcataa', result.sequence)
    end

    def test_invalid_arguments_raise
      assert_raise(ArgumentError) do
        Bio::CodonOptimizer.new(codon_table: 999)
      end
      assert_raise(ArgumentError) do
        Bio::CodonOptimizer.new(codon_table: 1,
                                forbidden_substrings: ['a?'])
      end
    end

    private

    def choices_for(sequence)
      codon_table = Bio::CodonTable[1]
      sequence.downcase.scan(/.{3}/).map do |codon|
        amino_acid = codon_table[codon]
        codon_table.revtrans(amino_acid).sort
      end
    end

    def brute_force(sequence, weights, forbidden, fixed_range,
                    gc_range = nil)
      fixed_positions = fixed_range.to_a
      codons = sequence.downcase.scan(/.{3}/)
      best = nil

      product(@choices).each do |candidate_codons|
        fixed_positions.each do |position|
          candidate_codons[position] = codons[position]
        end
        candidate = candidate_codons.join
        next if forbidden.any? { |pattern| candidate.include?(pattern) }
        next if gc_range && !gc_range.cover?(gc_count(candidate))

        score = candidate_codons.sum { |codon| weights.fetch(codon, 0) }
        if best && (score < best[:score] ||
                    (score == best[:score] &&
                     candidate >= best[:sequence]))
          next
        end

        best = { sequence: candidate, score: score }
      end

      best
    end

    def product(arrays)
      arrays.inject([[]]) do |combined, values|
        combined.flat_map do |prefix|
          values.map { |value| prefix + [value] }
        end
      end
    end

    def gc_count(sequence)
      sequence.count('g') + sequence.count('c')
    end

    def translate(sequence, table = Bio::CodonTable[1])
      Bio::Sequence::NA.new(sequence).translate(1, table).to_s
    end
  end
end
