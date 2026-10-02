# frozen_string_literal: true
#
# test/unit/bio/util/test_rna_secondary_structure.rb - Unit test for
# Bio::RNASecondaryStructure.
#
# Copyright::  Copyright (C) 2026 BioRuby Project
# License::    The Ruby License
#

# loading helper routine for testing bioruby
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

# libraries needed for the tests
require 'test/unit'
require 'bio/util/rna_secondary_structure'

module Bio
  class TestRNASecondaryStructureBasic < Test::Unit::TestCase
    def test_simple_hairpin_default_weights
      result = Bio::RNASecondaryStructure.new('GGGAAACCC', min_hairpin: 3).fold
      assert_equal(true, result.feasible?)
      assert_equal(3, result.score)
      assert_equal(1, result.count)
      assert_equal([[0, 8], [1, 7], [2, 6]], result.pairs)
      assert_equal('(((...)))', result.dot_bracket)
    end

    def test_class_method_fold
      result = Bio::RNASecondaryStructure.fold('GGGAAACCC', min_hairpin: 3)
      assert_equal(3, result.score)
      assert_equal('(((...)))', result.dot_bracket)
    end

    def test_dna_sequence_treated_as_rna
      result = Bio::RNASecondaryStructure.new('GGGTTTCCC', min_hairpin: 3).fold
      assert_equal(3, result.score)
      assert_equal('(((...)))', result.dot_bracket)
    end

    def test_empty_sequence
      result = Bio::RNASecondaryStructure.new('').fold
      assert_equal(true, result.feasible?)
      assert_equal(0, result.score)
      assert_equal(1, result.count)
      assert_equal([], result.pairs)
      assert_equal('', result.dot_bracket)
    end

    def test_multiple_optimal_structures_are_counted_once
      # two optimal structures: {(0,1),(2,3)} and {(0,3),(1,2)}
      result = Bio::RNASecondaryStructure.new('GCGC', min_hairpin: 0).fold
      assert_equal(2, result.score)
      assert_equal(2, result.count)
    end

    def test_custom_weight_block
      seq = 'GAAAC'
      weights = { 'GC' => 5 }
      result = Bio::RNASecondaryStructure.new(seq, min_hairpin: 1) do |i, j|
        weights[seq[i, 1] + seq[j, 1]]
      end.fold
      assert_equal(5, result.score)
      assert_equal(1, result.count)
      assert_equal([[0, 4]], result.pairs)
    end

    def test_weight_block_receives_zero_based_positions
      result = Bio::RNASecondaryStructure.new('ACGU', min_hairpin: 0) do |i, j|
        i == 0 && j == 3 ? 10 : nil
      end.fold
      assert_equal(10, result.score)
      assert_equal([[0, 3]], result.pairs)
    end
  end

  class TestRNASecondaryStructureConstraints < Test::Unit::TestCase
    def test_forbidden_pair_is_avoided
      result = Bio::RNASecondaryStructure.new('GGGAAAACCCC', min_hairpin: 3,
                                                            forbidden: [[2, 8]]).fold
      assert_equal(3, result.score)
      assert_equal([[0, 9], [1, 8], [2, 7]], result.pairs)
      assert_equal(false, result.pairs.include?([2, 8]))
    end

    def test_forced_pair_is_present
      result = Bio::RNASecondaryStructure.new('GGGAAAACCCC', min_hairpin: 3,
                                                            forced: [[0, 7]]).fold
      assert_equal(1, result.score)
      assert_equal([[0, 7]], result.pairs)
      assert_equal(1, result.count)
    end

    def test_min_hairpin_is_respected
      result = Bio::RNASecondaryStructure.new('GGGAAACCC', min_hairpin: 4).fold
      assert_equal(2, result.score)
      # five optimal structures: (0,7)+(1,6), (0,8)+(1,6), (0,8)+(1,7),
      # (0,8)+(2,7) and (1,8)+(2,7)
      assert_equal(5, result.count)
      assert_equal([[1, 8], [2, 7]], result.pairs)
      result.pairs.each do |i, j|
        assert(j - i - 1 >= 4, "hairpin too small for pair (#{i}, #{j})")
      end
    end

    def test_infeasible_when_forced_pair_has_no_weight
      result = Bio::RNASecondaryStructure.new('GGGAAACCC', min_hairpin: 0,
                                                            forced: [[0, 4]]).fold
      assert_equal(false, result.feasible?)
      assert_equal(0, result.count)
      assert_equal(-Float::INFINITY, result.score)
      assert_nil(result.pairs)
      assert_nil(result.dot_bracket)
    end

    def test_pairs_are_normalized
      folder = Bio::RNASecondaryStructure.new('GGGAAACCC', min_hairpin: 3,
                                                           forced: [[8, 0]],
                                                           forbidden: [[7, 1]])
      assert_equal([[0, 8]], folder.forced_pairs)
      assert_equal([[1, 7]], folder.forbidden_pairs)
    end
  end

  class TestRNASecondaryStructureConflicts < Test::Unit::TestCase
    def test_crossing_forced_pairs
      error = assert_raise(ArgumentError) do
        Bio::RNASecondaryStructure.new('AAAAAAAA', min_hairpin: 0,
                                                 forced: [[0, 5], [2, 7]])
      end
      assert_match(/forced pair \(0, 5\)/, error.message)
      assert_match(/forced pair \(2, 7\)/, error.message)
      assert_match(/cross/, error.message)
    end

    def test_forced_pairs_sharing_a_base
      error = assert_raise(ArgumentError) do
        Bio::RNASecondaryStructure.new('AAAAAAAAA', min_hairpin: 0,
                                                  forced: [[0, 5], [5, 8]])
      end
      assert_match(/forced pair \(0, 5\)/, error.message)
      assert_match(/forced pair \(5, 8\)/, error.message)
      assert_match(/position 5/, error.message)
    end

    def test_forced_pair_conflicts_with_forbidden_pair
      error = assert_raise(ArgumentError) do
        Bio::RNASecondaryStructure.new('AAAAAAAAA', min_hairpin: 0,
                                                  forced: [[0, 8]],
                                                  forbidden: [[0, 8]])
      end
      assert_match(/forced pair \(0, 8\)/, error.message)
      assert_match(/forbidden pair \(0, 8\)/, error.message)
    end

    def test_forced_pair_conflicts_with_min_hairpin
      error = assert_raise(ArgumentError) do
        Bio::RNASecondaryStructure.new('AAAAAAAAA', min_hairpin: 3,
                                                  forced: [[1, 3]])
      end
      assert_match(/forced pair \(1, 3\)/, error.message)
      assert_match(/minimum hairpin length 3/, error.message)
    end

    def test_invalid_pair_positions
      assert_raise(ArgumentError) do
        Bio::RNASecondaryStructure.new('AAAA', forced: [[0, 100]])
      end
      assert_raise(ArgumentError) do
        Bio::RNASecondaryStructure.new('AAAA', forced: [[-1, 2]])
      end
      assert_raise(ArgumentError) do
        Bio::RNASecondaryStructure.new('AAAA', forced: [[2, 2]])
      end
      assert_raise(ArgumentError) do
        Bio::RNASecondaryStructure.new('AAAA', forced: [1])
      end
    end

    def test_invalid_min_hairpin
      assert_raise(ArgumentError) do
        Bio::RNASecondaryStructure.new('AAAA', min_hairpin: -1)
      end
    end
  end

  class TestRNASecondaryStructureBruteForce < Test::Unit::TestCase
    # Enumerates every legal pair set (no shared positions, no crossing
    # pairs, hairpin and forbidden constraints satisfied, all forced
    # pairs covered) and returns [best score, number of optimal sets].
    # This is deliberately independent of the dynamic programming
    # recursion used by Bio::RNASecondaryStructure.
    def brute_force(length, min_hairpin, forced, forbidden, weight)
      candidates = []
      (0...length).each do |i|
        (i + min_hairpin + 1...length).each do |j|
          next if forbidden.include?([i, j])

          w = weight.call(i, j)
          next if w.nil? || w == false

          candidates << [i, j, w]
        end
      end
      best = -Float::INFINITY
      count = 0
      walk = nil
      walk = lambda do |idx, chosen, score|
        if idx == candidates.length
          covered = forced.all? { |f| chosen.any? { |c| c[0] == f[0] && c[1] == f[1] } }
          if covered
            if score > best
              best = score
              count = 1
            elsif score == best
              count += 1
            end
          end
          next
        end
        i, j, w = candidates[idx]
        walk.call(idx + 1, chosen + [[i, j]], score + w) if compatible?(chosen, i, j)
        walk.call(idx + 1, chosen, score)
      end
      walk.call(0, [], 0)
      [best, count]
    end

    # True when pair (i, j) can be added to the chosen pairs without
    # sharing a position and without crossing any of them.
    def compatible?(chosen, i, j)
      chosen.none? { |a, b| pair_conflict?(a, b, i, j) }
    end

    def pair_conflict?(a, b, i, j)
      return true if a == i || a == j || b == i || b == j

      (a < i && i < b && b < j) || (i < a && a < j && j < b)
    end

    def structure_score(pairs, weight)
      pairs.sum { |i, j| weight.call(i, j) }
    end

    def assert_valid_structure(result, length, min_hairpin, forced, forbidden, weight)
      pairs = result.pairs
      assert_equal(pairs.sort, pairs)
      pairs.each do |i, j|
        assert(j - i - 1 >= min_hairpin, "hairpin too small for pair (#{i}, #{j})")
        assert_equal(false, forbidden.include?([i, j]))
        assert_equal(false, weight.call(i, j).nil?)
      end
      pairs.combination(2) do |(a, b), (i, j)|
        assert_equal(false, pair_conflict?(a, b, i, j))
      end
      forced.each do |f|
        assert(pairs.include?(f), "forced pair #{f.inspect} missing")
      end
      assert_equal(result.score, structure_score(pairs, weight))
      assert_equal(length, result.dot_bracket.length)
    end

    def test_random_instances_against_brute_force
      rng = Random.new(20_261_002)
      200.times do |trial|
        length = rng.rand(0..12)
        seq = Array.new(length) { 'ACGU'[rng.rand(4)] }.join
        min_hairpin = rng.rand(0..3)
        weights = {}
        (0...length).each do |i|
          (i + 1...length).each do |j|
            weights[[i, j]] = rng.rand(1..3) if rng.rand < 0.5
          end
        end
        weight = ->(i, j) { weights[[i, j]] }
        forbidden = weights.keys.select { rng.rand < 0.1 }
        forced = pick_compatible_forced(rng, length, min_hairpin, forbidden)
        result = Bio::RNASecondaryStructure.new(seq, min_hairpin: min_hairpin,
                                                     forced: forced,
                                                     forbidden: forbidden,
                                                     &weight).fold
        best, count = brute_force(length, min_hairpin, forced, forbidden, weight)
        context = "seq=#{seq} min_hairpin=#{min_hairpin} forced=#{forced} " +                  "forbidden=#{forbidden} trial=#{trial}"
        assert_equal(best, result.score, context)
        assert_equal(count, result.count, context)
        if count.zero?
          assert_equal(false, result.feasible?, context)
          assert_nil(result.pairs, context)
        else
          assert_valid_structure(result, length, min_hairpin, forced,
                                 forbidden, weight)
        end
      end
    end

    # Picks up to two mutually compatible forced pairs (non-crossing, no
    # shared positions) that satisfy the hairpin and forbidden
    # constraints.  Pairs may have no weight, leading to infeasible
    # instances, which the brute force check also covers.
    def pick_compatible_forced(rng, length, min_hairpin, forbidden)
      candidates = []
      (0...length).each do |i|
        (i + min_hairpin + 1...length).each do |j|
          candidates << [i, j] unless forbidden.include?([i, j])
        end
      end
      candidates.shuffle!(random: rng)
      forced = []
      candidates.each do |i, j|
        break if forced.length >= 2
        next if forced.any? { |a, b| pair_conflict?(a, b, i, j) }

        forced << [i, j]
      end
      forced
    end
  end
end
