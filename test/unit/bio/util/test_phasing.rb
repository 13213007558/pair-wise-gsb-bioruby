# frozen_string_literal: true
#
# test/unit/bio/util/test_phasing.rb - Unit test for Bio::Phasing
#
# Copyright::  Copyright (C) 2026 BioRuby project
# License::    The Ruby License
#

# loading helper routine for testing bioruby
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

# libraries needed for the tests
require 'test/unit'
require 'bio/util/phasing'

module Bio
  class TestPhasingValidation < Test::Unit::TestCase
    def test_negative_weight_raises
      assert_raise(ArgumentError) do
        Bio::Phasing.new([{ 0 => 0 }], weights: [-1])
      end
    end

    def test_non_integer_weight_raises
      assert_raise(ArgumentError) do
        Bio::Phasing.new([{ 0 => 0 }], weights: [1.5])
      end
    end

    def test_weights_length_mismatch_raises
      assert_raise(ArgumentError) do
        Bio::Phasing.new([{ 0 => 0 }], weights: [1, 1])
      end
    end

    def test_invalid_allele_raises
      assert_raise(ArgumentError) { Bio::Phasing.new([{ 0 => 2 }]) }
    end

    def test_invalid_site_raises
      assert_raise(ArgumentError) { Bio::Phasing.new([{ '0' => 0 }]) }
    end

    def test_invalid_calls_type_raises
      assert_raise(ArgumentError) { Bio::Phasing.new([42]) }
    end

    def test_invalid_budget_raises
      assert_raise(ArgumentError) do
        Bio::Phasing.new([{ 0 => 0 }]).phase(max_nodes: -1)
      end
    end
  end

  class TestPhasing < Test::Unit::TestCase
    def test_perfect_fragments
      result = Bio::Phasing.new([
                                  { 0 => 0, 1 => 0, 2 => 1 },
                                  { 1 => 0, 2 => 1 },
                                  { 0 => 1, 1 => 1, 2 => 0 },
                                  { 1 => 1 }
                                ]).phase
      assert(result.optimal?)
      assert_equal(0, result.cost)
      assert_equal(0, result.lower_bound)
      assert_equal([[0, 0, 1], [1, 1, 0]], result.haplotypes)
      assert_equal([0, 0, 1, 1], result.assignment)
    end

    def test_error_correction_cost
      # The first fragment carries one sequencing error.
      result = Bio::Phasing.new([
                                  { 0 => 0, 1 => 1, 2 => 1 },
                                  { 0 => 0, 1 => 0 },
                                  { 1 => 1, 2 => 0 },
                                  { 0 => 1, 2 => 0 }
                                ]).phase
      assert(result.optimal?)
      assert_equal(1, result.cost)
      assert_equal(1, result.lower_bound)
    end

    def test_weights_change_optimum
      fragments = [
        { 0 => 0, 1 => 0 },
        { 0 => 1, 1 => 1 },
        { 0 => 0, 1 => 1 }
      ]
      light = Bio::Phasing.new(fragments, weights: [10, 10, 1]).phase
      assert_equal(1, light.cost)
      heavy = Bio::Phasing.new(fragments, weights: [10, 10, 100]).phase
      assert_equal(10, heavy.cost)
      assert_equal([[0, 1], [1, 0]], heavy.haplotypes)
    end

    def test_missing_calls
      result = Bio::Phasing.new([
                                  { 0 => 0, 2 => 1 },
                                  { 0 => 0, 1 => nil, 2 => 1 },
                                  { 0 => 1, 1 => 1, 2 => 0 }
                                ]).phase
      assert(result.optimal?)
      assert_equal(0, result.cost)
      assert_equal([0, 1, 2], result.sites)
    end

    def test_array_of_pairs_input
      result = Bio::Phasing.new([
                                  [[0, 0], [1, 0]],
                                  [[0, 1], [1, 1]]
                                ]).phase
      assert(result.optimal?)
      assert_equal(0, result.cost)
    end

    def test_disconnected_blocks_solved_separately
      result = Bio::Phasing.new([
                                  { 0 => 0, 1 => 1 },
                                  { 0 => 1, 1 => 0 },
                                  { 5 => 1, 6 => 1 },
                                  { 5 => 0, 6 => 0 }
                                ]).phase
      assert(result.optimal?)
      assert_equal(0, result.cost)
      assert_equal(2, result.blocks.length)
      assert_equal([0, 1], result.blocks[0].sites)
      assert_equal([5, 6], result.blocks[1].sites)
      assert_equal([0, 1, 5, 6], result.sites)
      assert_equal([0, 1, 0, 1], result.assignment)
    end

    def test_fragment_association_beats_column_voting
      # Per-column majority voting would phase the sites as 01/10,
      # but the fragments are perfectly explained by 00/11.
      result = Bio::Phasing.new([
                                  { 0 => 0, 1 => 0 },
                                  { 0 => 1, 1 => 1 },
                                  { 0 => 0 },
                                  { 1 => 1 }
                                ], weights: [1, 1, 2, 2]).phase
      assert(result.optimal?)
      assert_equal(0, result.cost)
      assert_equal([[0, 0], [1, 1]], result.haplotypes)
    end

    def test_haplotypes_are_complementary
      result = Bio::Phasing.new([
                                  { 0 => 0, 1 => 1, 3 => 0 },
                                  { 0 => 1, 2 => 0 },
                                  { 1 => 0, 2 => 1, 3 => 1 }
                                ]).phase
      h1, h2 = result.haplotypes
      h1.zip(h2).each { |a, b| assert_equal(1, a + b) }
    end

    def test_empty_input
      result = Bio::Phasing.new([]).phase
      assert(result.optimal?)
      assert_equal(0, result.cost)
      assert_equal([], result.sites)
      assert_equal([[], []], result.haplotypes)
      assert_equal([], result.assignment)
    end

    def test_fragment_without_calls
      result = Bio::Phasing.new([{}, { 0 => 0 }, { 0 => 1 }]).phase
      assert_equal(0, result.cost)
      assert_equal(0, result.assignment[0])
    end
  end

  class TestPhasingBudget < Test::Unit::TestCase
    FRAGMENTS = [
      { 0 => 0, 1 => 0 },
      { 0 => 1, 1 => 1 },
      { 0 => 0, 1 => 1 }
    ].freeze

    def cost_of(result, fragments, weights)
      position = {}
      result.sites.each_with_index { |site, i| position[site] = i }
      h1, h2 = result.haplotypes
      total = 0
      fragments.each_with_index do |calls, i|
        hap = result.assignment[i].zero? ? h1 : h2
        calls.each do |site, allele|
          total += weights[i] if hap[position[site]] != allele
        end
      end
      total
    end

    def test_budget_exhaustion_is_reported
      weights = [1, 1, 1]
      result = Bio::Phasing.new(FRAGMENTS).phase(max_nodes: 0)
      assert(!result.optimal?)
      assert(result.lower_bound <= result.cost)
      assert_equal(cost_of(result, FRAGMENTS, weights), result.cost)
    end

    def test_unbounded_search_proves_optimum
      weights = [1, 1, 1]
      proven = Bio::Phasing.new(FRAGMENTS).phase(max_nodes: nil)
      assert(proven.optimal?)
      assert_equal(proven.cost, proven.lower_bound)
      assert_equal(cost_of(proven, FRAGMENTS, weights), proven.cost)
      limited = Bio::Phasing.new(FRAGMENTS).phase(max_nodes: 0)
      assert(limited.cost >= proven.cost)
      assert(limited.lower_bound <= proven.cost)
    end

    def test_default_budget_solves_small_instances
      result = Bio::Phasing.new(FRAGMENTS).phase
      assert(result.optimal?)
    end
  end

  class TestPhasingExhaustive < Test::Unit::TestCase
    def test_optimal_count_unique
      result = Bio::Phasing.new([{ 0 => 0 }, { 0 => 1 }]).exhaustive
      assert_equal(0, result.cost)
      assert_equal(1, result.optimal_count)
    end

    def test_optimal_count_modulo_haplotype_swap
      # The zero-weight fragment may go to either side, giving exactly
      # 2 optimal solutions. Swapping the two haplotypes does not
      # double the count.
      result = Bio::Phasing.new(
        [{ 0 => 0, 1 => 0 }, { 0 => 1, 1 => 1 }, { 0 => 0, 1 => 1 }],
        weights: [1, 1, 0]
      ).exhaustive
      assert_equal(0, result.cost)
      assert_equal(2, result.optimal_count)
    end

    def test_exhaustive_raises_for_large_blocks
      fragments = (0..24).map { |i| { 0 => i % 2, 1 => (i + 1) % 2 } }
      assert_raise(ArgumentError) do
        Bio::Phasing.new(fragments).exhaustive(max_fragments: 20)
      end
    end

    def test_class_method_exhaustive
      result = Bio::Phasing.exhaustive([{ 0 => 0 }, { 0 => 1 }])
      assert_equal(0, result.cost)
    end

    def random_instance(rng)
      n_sites = rng.rand(2..5)
      fragments = []
      weights = []
      rng.rand(2..8).times do
        calls = {}
        n_sites.times { |s| calls[s] = rng.rand(2) if rng.rand < 0.7 }
        fragments << calls
        weights << rng.rand(0..3)
      end
      [fragments, weights]
    end

    def test_exact_search_matches_exhaustive_on_small_instances
      rng = Random.new(20_261_001)
      30.times do
        fragments, weights = random_instance(rng)
        phasing = Bio::Phasing.new(fragments, weights: weights)
        result = phasing.phase(max_nodes: nil)
        exhaustive = phasing.exhaustive
        assert(result.optimal?)
        assert_equal(exhaustive.cost, result.cost)
        assert_equal(result.cost, result.lower_bound)
      end
    end
  end
end
