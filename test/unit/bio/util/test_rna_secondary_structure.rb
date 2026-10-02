# frozen_string_literal: true

#
# test/unit/bio/util/test_rna_secondary_structure.rb -
#   Unit test for Bio::RNASecondaryStructure
#

# loading helper routine for testing bioruby
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

# libraries needed for the tests
require 'test/unit'
require 'bio/sequence'
require 'bio/util/rna_secondary_structure'

module Bio
  class TestRNASecondaryStructure < Test::Unit::TestCase
    def setup
      @klass = Bio::RNASecondaryStructure
    end

    def test_fold_default_scoring
      result = @klass.new('GGGAAACCC', min_loop: 3).fold
      assert_equal(3, result.score)
      assert_equal([[0, 8], [1, 7], [2, 6]], result.pairs)
      assert_equal('(((...)))', result.dot_bracket)
      assert_equal(1, result.count)
      assert_predicate(result, :feasible?)
    end

    def test_class_method_fold
      result = @klass.fold('GGGAAACCC', min_loop: 3)
      assert_equal(3, result.score)
      assert_equal('(((...)))', result.dot_bracket)
    end

    def test_accepts_sequence_object
      seq = Bio::Sequence::NA.new('gggaaaccc')
      result = @klass.new(seq, min_loop: 3).fold
      assert_equal(3, result.score)
    end

    def test_counts_all_optimal_structures
      # Every pair allowed, min_loop 0: the optima are
      # {(0,1),(2,3)} and {(0,3),(1,2)}.
      allow_all = ->(_i, _j, _s) { true }
      result = @klass.new('AAAA', min_loop: 0,
                                 pair_allowed: allow_all).fold
      assert_equal(2, result.score)
      assert_equal(2, result.count)
    end

    def test_empty_structure_counted_once
      # No pair can form: the empty structure is the single optimum.
      result = @klass.new('AAAAAA', min_loop: 3).fold
      assert_equal(0, result.score)
      assert_equal([], result.pairs)
      assert_equal('......', result.dot_bracket)
      assert_equal(1, result.count)
    end

    def test_min_loop_is_respected
      allow_all = ->(_i, _j, _s) { true }
      # With min_loop 2 the pair (0, 2) is too close; only (0, 3) and
      # (1, 3)-style pairs remain. Sequence length 4, min_loop 2:
      # allowed pairs need j - i > 2, so only (0, 3) exists.
      result = @klass.new('AAAA', min_loop: 2,
                                 pair_allowed: allow_all).fold
      assert_equal(1, result.score)
      assert_equal([[0, 3]], result.pairs)
      assert_equal(1, result.count)
    end

    def test_forced_pair_is_present
      # Forcing (0, 5) forbids the better nested chain inside.
      result = @klass.new('GAAAAC', min_loop: 3, forced: [[0, 5]]).fold
      assert_equal(1, result.score)
      assert_equal([[0, 5]], result.pairs)
      assert_equal('(....)', result.dot_bracket)
      assert_equal(1, result.count)
    end

    def test_forbidden_pair_is_avoided
      result = @klass.new('GGGAAACCC', min_loop: 3,
                                       forbidden: [[0, 8]]).fold
      assert_equal(2, result.score)
      assert_equal([[1, 7], [2, 6]], result.pairs)
      assert(!result.pairs.include?([0, 8]))
    end

    def test_custom_score_function
      # GC pairs score 3, AU pairs score 1.
      score = lambda do |i, j, s|
        pair = [s[i].upcase, s[j].upcase].sort.join
        pair == 'CG' ? 3 : 1
      end
      # GGGAAACCC: three GC pairs = 9; an AU-only alternative scores
      # less.
      result = @klass.new('GGGAAACCC', min_loop: 3, score: score).fold
      assert_equal(9, result.score)
      assert_equal(1, result.count)
    end

    def test_each_structure_enumerates_uniquely
      allow_all = ->(_i, _j, _s) { true }
      folder = @klass.new('AAAA', min_loop: 0, pair_allowed: allow_all)
      structures = folder.each_structure.to_a
      # {}, 6 single pairs, {(0,1),(2,3)}, {(0,3),(1,2)}
      assert_equal(9, structures.size)
      assert_equal(structures.uniq.size, structures.size)
      assert_include(structures, [])
      assert_include(structures, [[0, 1], [2, 3]])
      assert_include(structures, [[0, 3], [1, 2]])
    end

    def test_each_structure_respects_forced_pairs
      folder = @klass.new('GUAAC', min_loop: 3, forced: [[0, 4]])
      structures = folder.each_structure.to_a
      assert_equal(1, structures.size)
      assert_equal([[0, 4]], structures.first)
    end

    def test_result_to_s
      result = @klass.new('GGGAAACCC', min_loop: 3).fold
      assert_equal('(((...)))', result.to_s)
    end

    def test_empty_sequence
      result = @klass.new('', min_loop: 3).fold
      assert_equal(0, result.score)
      assert_equal([], result.pairs)
      assert_equal('', result.dot_bracket)
      assert_equal(1, result.count)
    end

    def test_conflicting_forced_pairs_crossing
      error = assert_raise(@klass::ConstraintConflictError) do
        @klass.new('G' * 12, min_loop: 0,
                            pair_allowed: ->(_i, _j, _s) { true },
                            forced: [[1, 8], [4, 10]])
      end
      assert_match(/forced pair \(1, 8\)/, error.message)
      assert_match(/forced pair \(4, 10\)/, error.message)
      assert_match(/cross/, error.message)
      assert_equal(2, error.constraints.size)
    end

    def test_conflicting_forced_pairs_shared_position
      error = assert_raise(@klass::ConstraintConflictError) do
        @klass.new('G' * 12, min_loop: 0,
                            pair_allowed: ->(_i, _j, _s) { true },
                            forced: [[1, 8], [8, 10]])
      end
      assert_match(/forced pair \(1, 8\)/, error.message)
      assert_match(/forced pair \(8, 10\)/, error.message)
      assert_match(/8/, error.message)
    end

    def test_conflicting_forced_and_forbidden
      error = assert_raise(@klass::ConstraintConflictError) do
        @klass.new('GAAAAC', min_loop: 3,
                            forced: [[0, 5]], forbidden: [[0, 5]])
      end
      assert_match(/forced pair \(0, 5\)/, error.message)
      assert_match(/forbidden pair \(0, 5\)/, error.message)
    end

    def test_conflicting_forced_pair_and_min_loop
      error = assert_raise(@klass::ConstraintConflictError) do
        @klass.new('GAC', min_loop: 3, forced: [[0, 2]])
      end
      assert_match(/forced pair \(0, 2\)/, error.message)
      assert_match(/hairpin loop/, error.message)
    end

    def test_conflicting_forced_pair_and_pairing_rules
      error = assert_raise(@klass::ConstraintConflictError) do
        @klass.new('AAAAAA', min_loop: 3, forced: [[0, 5]])
      end
      assert_match(/forced pair \(0, 5\)/, error.message)
      assert_match(/base-pairing rules/, error.message)
    end

    def test_out_of_range_pair_raises_argument_error
      assert_raise(ArgumentError) do
        @klass.new('GAAAAC', forced: [[0, 6]])
      end
      assert_raise(ArgumentError) do
        @klass.new('GAAAAC', forbidden: [[-1, 3]])
      end
    end

    def test_negative_min_loop_raises_argument_error
      assert_raise(ArgumentError) do
        @klass.new('GAAAAC', min_loop: -1)
      end
    end

    # --- cross-validation against exhaustive enumeration ------------

    def brute_force_optimum(folder, score_fn)
      best = nil
      count = 0
      best_structures = []
      folder.each_structure do |pairs|
        value = pairs.sum { |i, j| score_fn.call(i, j, folder.sequence) }
        if best.nil? || value > best
          best = value
          count = 1
          best_structures = [pairs]
        elsif value == best
          count += 1
          best_structures << pairs
        end
      end
      [best, count, best_structures]
    end

    def random_instance(rng, index)
      n = rng.rand(2..10)
      seq = Array.new(n) { 'ACGU'[rng.rand(4)] }.join
      min_loop = rng.rand(0..3)
      allowed = {}
      pair_allowed = lambda do |i, j, _s|
        key = [i, j]
        allowed[key] = rng.rand < 0.6 if allowed[key].nil?
        allowed[key]
      end
      scores = {}
      score = lambda do |i, j, _s|
        scores[[i, j]] ||= rng.rand(-2..3)
      end
      forced = pick_forced_pairs(rng, n, min_loop, pair_allowed)
      forbidden = []
      n.times do |a|
        n.times do |b|
          next unless a < b && b - a > min_loop
          next if forced.include?([a, b])
          forbidden << [a, b] if rng.rand < 0.15
        end
      end
      folder = @klass.new(seq, min_loop: min_loop, forced: forced,
                               forbidden: forbidden,
                               pair_allowed: pair_allowed, score: score)
      label = "instance #{index}: seq=#{seq} min_loop=#{min_loop} " \
              "forced=#{forced} forbidden=#{forbidden}"
      [folder, score, label]
    end

    def pick_forced_pairs(rng, n, min_loop, pair_allowed)
      candidates = []
      n.times do |a|
        n.times do |b|
          next unless a < b && b - a > min_loop
          candidates << [a, b] if pair_allowed.call(a, b, nil)
        end
      end
      candidates.shuffle!(random: rng)
      chosen = []
      candidates.each do |(a, b)|
        next if chosen.any? { |(c, d)| [c, d].intersect?([a, b]) }
        next if chosen.any? { |(c, d)| a < c && c < b && b < d }
        next if chosen.any? { |(c, d)| c < a && a < d && d < b }
        break if chosen.size >= 2

        chosen << [a, b]
      end
      chosen.sort
    end

    def test_matches_exhaustive_enumeration
      rng = Random.new(20_260_102)
      100.times do |index|
        folder, score_fn, label = random_instance(rng, index)
        result = folder.fold
        best, count, best_structures = brute_force_optimum(folder, score_fn)
        assert_equal(best, result.score, "score mismatch: #{label}")
        assert_equal(count, result.count, "count mismatch: #{label}")
        assert_include(best_structures, result.pairs,
                       "returned structure is not optimal: #{label}")
      end
    end

    def test_enumeration_has_no_duplicates
      rng = Random.new(917)
      20.times do |index|
        folder, _score_fn, label = random_instance(rng, index)
        structures = folder.each_structure.to_a
        assert_equal(structures.uniq.size, structures.size,
                     "duplicate structures enumerated: #{label}")
      end
    end
  end
end
