#
# test/unit/bio/util/test_haplotype_phaser.rb - Unit test for Bio::HaplotypePhaser.
#
# License:: The Ruby License
#

# loading helper routine for testing bioruby
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

# libraries needed for the tests
require 'test/unit'
require 'bio/util/haplotype_phaser'

module Bio
  class TestHaplotypePhaser < Test::Unit::TestCase
    def test_phase_consistent_fragments
      fragments = [
        { 0 => 0, 1 => 0 },
        { 0 => 1, 1 => 1 },
        { 1 => 0, 2 => 1 },
        { 2 => 0 },
      ]
      result = Bio::HaplotypePhaser.new(fragments).phase
      assert_equal(true, result.optimal?)
      assert_equal(0, result.cost)
      assert_equal(0, result.lower_bound)
      assert_equal(0, result.upper_bound)
      assert_equal([0, 1, 2], result.sites)
      assert_consistent(result, fragments)
    end

    def test_array_input_with_missing_calls
      fragments = [
        [0, 0, nil],
        [1, 1, nil],
        [nil, 0, 1],
      ]
      result = Bio::HaplotypePhaser.new(fragments).phase
      assert_equal(0, result.cost)
      assert_equal([0, 1, 2], result.sites)
      assert_equal(true, result.optimal?)
    end

    def test_known_optimum
      # Hand-computed: the last fragment conflicts with both haplotypes
      # on exactly one site, so the optimum is 1.
      fragments = [
        { 0 => 0, 1 => 0 },
        { 1 => 0, 2 => 1 },
        { 0 => 1, 1 => 1 },
        { 1 => 1, 2 => 0 },
        { 0 => 0, 2 => 0 },
      ]
      result = Bio::HaplotypePhaser.new(fragments).phase
      assert_equal(1, result.cost)
      assert_equal(true, result.optimal?)
      assert_equal(1, Bio::HaplotypePhaser.exhaustive(fragments).cost)
    end

    def test_weights_change_optimum
      # Three pairwise conflicting fragments: the optimum is the
      # cheapest column minority, which depends on the weights.
      fragments = [
        { 0 => 0, 1 => 0 },
        { 0 => 1, 1 => 1 },
        { 0 => 0, 1 => 1 },
      ]
      assert_equal(1, Bio::HaplotypePhaser.new(fragments).phase.cost)
      weighted = Bio::HaplotypePhaser.new(fragments,
                                          weights: [2, 3, 4])
      assert_equal(2, weighted.phase.cost)
      assert_equal(2, Bio::HaplotypePhaser.exhaustive(
        fragments, weights: [2, 3, 4]
      ).cost)
    end

    def test_disconnected_blocks_solved_separately
      fragments = [
        { 0 => 0, 1 => 1 },
        { 0 => 1, 1 => 0 },
        { 5 => 1 },
        { 5 => 1, 6 => 0 },
      ]
      result = Bio::HaplotypePhaser.new(fragments).phase
      assert_equal(2, result.blocks.size)
      assert_equal([[0, 1], [5, 6]], result.blocks.map(&:sites))
      assert_equal([0, 1, 5, 6], result.sites)
      assert_equal(0, result.cost)
      assert_equal(true, result.optimal?)
      assert_consistent(result, fragments)
    end

    def test_budget_limits_search
      fragments = random_fragments(Random.new(42), 10, 6, 0.2)
      limited = Bio::HaplotypePhaser.new(fragments, budget: 1).phase
      assert_equal(false, limited.optimal?)
      assert_operator(limited.lower_bound, :<=, limited.cost)
      assert_equal(limited.cost, limited.upper_bound)
      full = Bio::HaplotypePhaser.new(fragments).phase
      assert_equal(true, full.optimal?)
      assert_operator(limited.lower_bound, :<=, full.cost)
      assert_operator(limited.cost, :>=, full.cost)
    end

    def test_exhaustive_matches_branch_and_bound
      rng = Random.new(1234)
      25.times do
        fragments = random_fragments(rng, rng.rand(2..8), rng.rand(2..5), 0.3)
        weights = Array.new(fragments.size) { rng.rand(1..3) }
        result = Bio::HaplotypePhaser.new(fragments, weights: weights).phase
        exhaustive = Bio::HaplotypePhaser.exhaustive(fragments,
                                                     weights: weights)
        assert_equal(true, result.optimal?)
        assert_equal(exhaustive.cost, result.cost)
        assert_equal(result.cost, result.lower_bound)
      end
    end

    def test_exhaustive_counts_swap_as_one_solution
      # Two conflicting fragments have a unique optimal assignment once
      # swapping the two haplotypes is factored out.
      conflict = Bio::HaplotypePhaser.exhaustive([{ 0 => 0 }, { 0 => 1 }])
      assert_equal(0, conflict.cost)
      assert_equal(1, conflict.solutions)
      pairs = Bio::HaplotypePhaser.exhaustive(
        [{ 0 => 0 }, { 0 => 0 }, { 0 => 1 }, { 0 => 1 }]
      )
      assert_equal(0, pairs.cost)
      assert_equal(1, pairs.solutions)
    end

    def test_exhaustive_counts_distinct_optima
      # Two fragments covering disjoint sites: both assignments are
      # optimal, and swapping haplotypes does not add new solutions.
      result = Bio::HaplotypePhaser.exhaustive([{ 0 => 0 }, { 1 => 1 }])
      assert_equal(0, result.cost)
      assert_equal(2, result.solutions)
    end

    def test_exhaustive_limit
      fragments = Array.new(25) { { 0 => 0 } }
      assert_raise(ArgumentError) do
        Bio::HaplotypePhaser.exhaustive(fragments)
      end
    end

    def test_empty_input
      result = Bio::HaplotypePhaser.new([]).phase
      assert_equal(0, result.cost)
      assert_equal(true, result.optimal?)
      assert_equal([], result.assignments)
      exhaustive = Bio::HaplotypePhaser.exhaustive([])
      assert_equal(0, exhaustive.cost)
      assert_equal(1, exhaustive.solutions)
    end

    def test_validation
      assert_raise(ArgumentError) { Bio::HaplotypePhaser.new([{ 0 => 2 }]) }
      assert_raise(ArgumentError) { Bio::HaplotypePhaser.new([{ -1 => 0 }]) }
      assert_raise(ArgumentError) { Bio::HaplotypePhaser.new(['x']) }
      assert_raise(ArgumentError) do
        Bio::HaplotypePhaser.new([{ 0 => 0 }], weights: [-1])
      end
      assert_raise(ArgumentError) do
        Bio::HaplotypePhaser.new([{ 0 => 0 }], weights: [1, 2])
      end
      assert_raise(ArgumentError) do
        Bio::HaplotypePhaser.new([{ 0 => 0 }], budget: -1)
      end
    end

    private

    def random_fragments(rng, count, nsites, missing_prob)
      Array.new(count) do
        frag = {}
        nsites.times do |site|
          frag[site] = rng.rand(2) if rng.rand >= missing_prob
        end
        frag
      end
    end

    # With cost 0 every fragment must match its assigned haplotype on
    # all of its called sites.
    def assert_consistent(result, fragments)
      fragments.each_with_index do |frag, i|
        hap = result.haplotypes[result.assignments[i]]
        frag.each do |site, allele|
          pos = result.sites.index(site)
          assert_equal(allele, hap[pos],
                       "fragment #{i} mismatches at site #{site}")
        end
      end
    end
  end
end
