# frozen_string_literal: true
#
# test/unit/bio/util/test_rna_fold.rb - Unit test for Bio::RNAFold.
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
require 'bio/util/rna_fold'

module Bio
  class TestRNAFold < Test::Unit::TestCase
    # Independent reference: enumerates every non-crossing matching of
    # 0...n whose pairs span at least min_hairpin + 1 positions and are
    # allowed by the energy hash, and computes ensemble averages directly.
    class BruteForce
      attr_reader :partition_function, :pair_probabilities, :unpaired_probabilities

      def initialize(n, energies, min_hairpin,
                     forced_pairs: [], forbidden_pairs: [], forced_unpaired: [])
        @n = n
        @energies = energies
        @min_hairpin = min_hairpin
        structures = enumerate(0, n - 1).select do |pairs|
          forced_pairs.all? { |p| pairs.include?(p) } &&
            forbidden_pairs.none? { |p| pairs.include?(p) } &&
            forced_unpaired.all? { |pos| pairs.none? { |pr| pr.include?(pos) } }
        end
        @partition_function = structures.sum { |pairs| weight(pairs) }
        @pair_probabilities = Hash.new(0.0)
        @unpaired_probabilities = Array.new(n, 0.0)
        structures.each do |pairs|
          w = weight(pairs)
          pairs.each { |p| @pair_probabilities[p] += w }
          (0...n).each do |pos|
            @unpaired_probabilities[pos] += w if pairs.none? { |pr| pr.include?(pos) }
          end
        end
        z = @partition_function
        @pair_probabilities.each_key { |p| @pair_probabilities[p] /= z }
        @unpaired_probabilities.map! { |x| x / z }
      end

      private

      def weight(pairs)
        pairs.inject(1.0) { |prod, (i, j)| prod * Math.exp(-@energies[[i, j]]) }
      end

      def enumerate(lo, hi)
        return [[]] if lo > hi
        structures = enumerate(lo + 1, hi)
        (lo + @min_hairpin + 1..hi).each do |k|
          next unless @energies[[lo, k]]
          enumerate(lo + 1, k - 1).each do |inside|
            enumerate(k + 1, hi).each do |right|
              structures << ([[lo, k]] + inside + right)
            end
          end
        end
        structures
      end
    end

    def random_energies(n, min_hairpin, rng, energy_range)
      energies = {}
      (0...n).each do |i|
        (i + min_hairpin + 1...n).each do |j|
          next unless rng.rand < 0.6
          energies[[i, j]] = energy_range.begin +
                             rng.rand * (energy_range.end - energy_range.begin)
        end
      end
      energies
    end

    def assert_matches_brute_force(fold, energies, min_hairpin, constraints = {})
      n = fold.length
      ref = BruteForce.new(n, energies, min_hairpin, **constraints)
      assert_in_delta(Math.log(ref.partition_function),
                      fold.log_partition_function, 1e-8,
                      'partition function mismatch')
      (0...n).each do |i|
        (i + 1...n).each do |j|
          key = [i, j]
          next unless energies[key]
          assert_in_delta(ref.pair_probabilities[key],
                          fold.pair_probability(i, j), 1e-8,
                          "pair probability (#{i}, #{j}) mismatch")
        end
        assert_in_delta(ref.unpaired_probabilities[i],
                        fold.unpaired_probability(i), 1e-8,
                        "unpaired probability #{i} mismatch")
        paired = (0...n).sum do |j|
          next 0.0 if i == j
          fold.pair_probability([i, j].min, [i, j].max)
        end
        assert_in_delta(1.0, paired + fold.unpaired_probability(i), 1e-8,
                        "probabilities for base #{i} do not sum to 1")
      end
    end

    def test_random_sequences_against_brute_force
      rng = Random.new(1234)
      10.times do
        n = 8 + rng.rand(4)
        min_hairpin = [0, 1, 3][rng.rand(3)]
        energies = random_energies(n, min_hairpin, rng, -4.0..2.0)
        fold = Bio::RNAFold.new('A' * n, pair_energy: energies,
                                         min_hairpin: min_hairpin)
        assert_matches_brute_force(fold, energies, min_hairpin)
      end
    end

    def test_constraints_against_brute_force
      rng = Random.new(42)
      n = 10
      min_hairpin = 1
      energies = random_energies(n, min_hairpin, rng, -3.0..1.0)
      constraints = {
        forced_pairs: [[1, 8]],
        forbidden_pairs: energies.keys.select { |i, _j| i == 0 },
        forced_unpaired: [4]
      }
      constraints[:forced_pairs].each { |p| energies[p] ||= -2.0 }
      fold = Bio::RNAFold.new('A' * n, pair_energy: energies,
                                       min_hairpin: min_hairpin, **constraints)
      assert_matches_brute_force(fold, energies, min_hairpin, constraints)
    end

    def test_tiny_example_partition_function
      # Pairs (0,2), (0,3), (1,3) with min_hairpin 1; no two can coexist.
      energies = { [0, 2] => -1.0, [0, 3] => -2.0, [1, 3] => -0.5 }
      fold = Bio::RNAFold.new('AAAA', pair_energy: energies, min_hairpin: 1)
      expected = 1.0 + Math.exp(1.0) + Math.exp(2.0) + Math.exp(0.5)
      assert_in_delta(expected, fold.partition_function, 1e-10)
    end

    def test_forced_pair_has_probability_one
      energies = { [0, 4] => -2.0, [1, 3] => -1.0, [0, 3] => -1.0 }
      fold = Bio::RNAFold.new('AAAAA', pair_energy: energies, min_hairpin: 1,
                                        forced_pairs: [[0, 4]])
      assert_in_delta(1.0, fold.pair_probability(0, 4), 1e-12)
      assert_in_delta(0.0, fold.unpaired_probability(0), 1e-12)
      assert_in_delta(0.0, fold.unpaired_probability(4), 1e-12)
    end

    def test_forbidden_pair_has_probability_zero
      energies = { [0, 3] => -2.0, [1, 3] => -1.0 }
      fold = Bio::RNAFold.new('AAAA', pair_energy: energies, min_hairpin: 1,
                                       forbidden_pairs: [[0, 3]])
      assert_equal(0.0, fold.pair_probability(0, 3))
      assert(fold.pair_probability(1, 3) > 0.0)
    end

    def test_forced_unpaired_base
      energies = { [0, 3] => -2.0, [1, 3] => -1.0 }
      fold = Bio::RNAFold.new('AAAA', pair_energy: energies, min_hairpin: 1,
                                       forced_unpaired: [3])
      assert_equal(1.0, fold.unpaired_probability(3))
      assert_equal(0.0, fold.pair_probability(0, 3))
      assert_equal(0.0, fold.pair_probability(1, 3))
    end

    def test_conflicting_constraints_are_rejected
      energies = { [0, 4] => -2.0, [1, 3] => -1.0, [2, 6] => -1.0,
                   [0, 6] => -1.0, [1, 4] => -1.0 }
      # forced pair that is also forbidden
      assert_raise(Bio::RNAFold::ConstraintError) do
        Bio::RNAFold.new('AAAAAAA', pair_energy: energies, min_hairpin: 1,
                                    forced_pairs: [[0, 4]], forbidden_pairs: [[0, 4]])
      end
      # forced pair whose base must be unpaired
      assert_raise(Bio::RNAFold::ConstraintError) do
        Bio::RNAFold.new('AAAAAAA', pair_energy: energies, min_hairpin: 1,
                                    forced_pairs: [[0, 4]], forced_unpaired: [4])
      end
      # two forced pairs sharing a base
      assert_raise(Bio::RNAFold::ConstraintError) do
        Bio::RNAFold.new('AAAAAAA', pair_energy: energies, min_hairpin: 1,
                                    forced_pairs: [[0, 4], [4, 6]])
      end
      # crossing forced pairs (a pseudoknot), in either order
      assert_raise(Bio::RNAFold::ConstraintError) do
        Bio::RNAFold.new('AAAAAAA', pair_energy: energies, min_hairpin: 1,
                                    forced_pairs: [[0, 4], [1, 6]])
      end
      assert_raise(Bio::RNAFold::ConstraintError) do
        Bio::RNAFold.new('AAAAAAA', pair_energy: energies, min_hairpin: 1,
                                    forced_pairs: [[1, 6], [0, 4]])
      end
      # forced pair violating min_hairpin
      assert_raise(Bio::RNAFold::ConstraintError) do
        Bio::RNAFold.new('AAAAAAA', pair_energy: energies, min_hairpin: 3,
                                    forced_pairs: [[1, 3]])
      end
      # forced pair that the energy model does not allow
      assert_raise(Bio::RNAFold::ConstraintError) do
        Bio::RNAFold.new('AAAAAAA', pair_energy: energies, min_hairpin: 1,
                                    forced_pairs: [[3, 6]])
      end
    end

    def test_nested_forced_pairs_are_allowed
      energies = { [0, 6] => -2.0, [2, 4] => -1.0 }
      fold = Bio::RNAFold.new('AAAAAAA', pair_energy: energies, min_hairpin: 1,
                                          forced_pairs: [[0, 6], [2, 4]])
      assert_in_delta(1.0, fold.pair_probability(0, 6), 1e-12)
      assert_in_delta(1.0, fold.pair_probability(2, 4), 1e-12)
    end

    def test_numerical_stability_with_extreme_energies
      n = 12
      energies = {}
      (0...n).each do |i|
        (i + 4...n).each do |j|
          energies[[i, j]] = (i + j).even? ? -600.0 : 300.0
        end
      end
      fold = Bio::RNAFold.new('A' * n, pair_energy: energies, min_hairpin: 3)
      assert(fold.log_partition_function.finite?)
      (0...n).each do |i|
        paired = (0...n).sum do |j|
          next 0.0 if i == j
          fold.pair_probability([i, j].min, [i, j].max)
        end
        unpaired = fold.unpaired_probability(i)
        assert(unpaired >= 0.0 && unpaired <= 1.0)
        assert_in_delta(1.0, paired + unpaired, 1e-8)
      end
    end

    def test_no_allowed_pairs
      fold = Bio::RNAFold.new('AAAAAA', pair_energy: {}, min_hairpin: 1)
      assert_in_delta(0.0, fold.log_partition_function, 1e-12)
      assert_equal([1.0] * 6, fold.unpaired_probabilities)
      assert_equal({}, fold.pair_probabilities)
    end

    def test_pair_energy_forms
      hash = { [0, 3] => -2.0 }
      matrix = Array.new(4) { Array.new(4) }
      matrix[0][3] = -2.0
      callable = ->(i, j) { [i, j] == [0, 3] ? -2.0 : nil }
      [hash, matrix, callable].each do |source|
        fold = Bio::RNAFold.new('AAAA', pair_energy: source, min_hairpin: 1)
        assert_equal(-2.0, fold.pair_energy(0, 3))
        assert_nil(fold.pair_energy(0, 2))
        assert_in_delta(Math.exp(2.0), fold.pair_weight(0, 3), 1e-10)
      end
    end

    def test_accepts_sequence_na
      require 'bio/sequence'
      seq = Bio::Sequence::NA.new('GGGAAACCC')
      fold = Bio::RNAFold.new(seq, pair_energy: { [0, 8] => -2.0 },
                                   min_hairpin: 3)
      assert_equal(9, fold.length)
      assert(fold.pair_probability(0, 8) > 0.0)
    end

    def test_invalid_arguments
      assert_raise(ArgumentError) do
        Bio::RNAFold.new('AAAA', pair_energy: {}, min_hairpin: -1)
      end
      assert_raise(ArgumentError) do
        Bio::RNAFold.new('AAAA', pair_energy: 'not an energy source')
      end
      assert_raise(ArgumentError) do
        Bio::RNAFold.new('AAAA', pair_energy: { [0, 3] => -1.0 },
                                  forced_pairs: [[0, 10]])
      end
      assert_raise(ArgumentError) do
        Bio::RNAFold.new('AAAA', pair_energy: { [0, 3] => Float::NAN },
                                  min_hairpin: 1)
      end
    end
  end
end
