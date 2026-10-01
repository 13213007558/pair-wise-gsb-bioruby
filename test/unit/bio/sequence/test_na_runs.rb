# frozen_string_literal: true

#
# = test/unit/bio/sequence/test_na_runs.rb - Unit test for tandem repeat runs
#
# Copyright::   Copyright (C) 2026 BioRuby Project
# License::     The Ruby License
#

# loading helper routine for testing bioruby
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 3,
                            'bioruby_test_helper.rb')).cleanpath.to_s

# libraries needed for the tests
require 'test/unit'
require 'bio/sequence'
require 'bio/sequence/na'

module Bio
  class TestSequenceNARuns < Test::Unit::TestCase
    def setup
      @run_class = Bio::Sequence::NA::Run
    end

    def run_array(sequence)
      Bio::Sequence::NA.new(sequence).runs.map(&:to_a)
    end

    def test_empty
      assert_equal([], run_array(''))
    end

    def test_no_repetition
      assert_equal([], run_array('acgt'))
    end

    def test_homopolymer
      assert_equal([[0, 4, 1]], run_array('aaaa'))
      runs = run_array('caaaag')
      assert_equal([[1, 5, 1]], runs)
    end

    def test_tandem_repeat
      assert_equal([[0, 6, 3]], run_array('atgatg'))
      assert_equal([[0, 8, 3]], run_array('atgatgat'))
    end

    def test_alternating
      # One run of period 2 for a maximal alternating segment.
      assert_equal([[0, 5, 2]], run_array('ababa'.tr('b', 't')))
    end

    def test_overlapping_runs
      # [0, 5) has period 2 and [4, 7) has period 1 and they overlap.
      assert_equal([[0, 5, 2], [4, 7, 1]], run_array('atataaa'))
    end

    def test_ambiguity_is_a_boundary
      # The ambiguous 'n' separates the two homopolymer runs.
      assert_equal([[0, 2, 1], [3, 5, 1]], run_array('aanaa'))
      # Every IUPAC ambiguity code is a boundary.
      'rywskmbdhvn'.each_char do |base|
        assert_equal([[0, 2, 1], [3, 5, 1]], run_array("aa#{base}aa"),
                     "#{base} should be a boundary")
      end
    end

    def test_u_is_concrete
      assert_equal([[0, 4, 2]], run_array('auau'))
    end

    def test_run_attributes
      run = @run_class.new(2, 9, 3)
      assert_equal(2, run.start_position)
      assert_equal(9, run.end_position)
      assert_equal(3, run.period)
      assert_equal(7, run.length)
      assert_equal(7, run.size)
      assert_in_delta(7.0 / 3, run.exponent, 1e-9)
      assert_equal([2, 9, 3], run.to_a)
      assert_equal(run, @run_class.new(2, 9, 3))
      assert_not_equal(run, @run_class.new(2, 8, 3))
      assert_equal('#<Bio::Sequence::NA::Run: [2, 9, 3]>',
                   run.inspect)
    end

    def test_sorted_by_position
      runs = run_array('aattaataattaa')
      starts = runs.map { |run| run[0] }
      assert_equal(starts.sort, starts)
    end

    # Exhaustive reference implementation: enumerates every interval,
    # determines its smallest period and checks maximality directly.
    def reference_runs(sequence)
      bytes = sequence.bytes
      length = bytes.length
      result = []
      (0...length).each do |start|
        (start + 2..length).each do |finish|
          period = smallest_period(bytes[start...finish])
          next if finish - start < 2 * period
          next if start > 0 &&
                  periodic?(bytes[(start - 1)...finish], period)
          next if finish < length &&
                  periodic?(bytes[start...(finish + 1)], period)
          result << [start, finish, period]
        end
      end
      result.sort
    end

    def smallest_period(bytes)
      z = z_function(bytes)
      length = bytes.length
      (1...length).each do |period|
        return period if z[period] >= length - period
      end
      length
    end

    def periodic?(bytes, period)
      (period...bytes.length).all? do |index|
        bytes[index] == bytes[index - period]
      end
    end

    def z_function(data)
      length = data.length
      z = Array.new(length, 0)
      left = 0
      right = 0
      (1...length).each do |index|
        z[index] = [right - index, z[index - left]].min if index < right
        z[index] += 1 while index + z[index] < length &&
                           data[z[index]] == data[index + z[index]]
        if index + z[index] > right
          left = index
          right = index + z[index]
        end
      end
      z
    end

    def assert_matches_reference(sequence)
      expected = reference_runs(sequence)
      actual = run_array(sequence)
      assert_equal(expected, actual,
                   "runs mismatch for #{sequence.inspect}")
    end

    def test_exhaustive_binary_alphabet
      # All binary strings up to length 7.
      1.upto(7) do |length|
        0.upto((1 << length) - 1) do |mask|
          sequence = length.times.map do |index|
            mask[index].zero? ? 'a' : 't'
          end.join
          assert_matches_reference(sequence)
        end
      end
    end

    def test_exhaustive_ternary_alphabet
      # All strings over three bases up to length 5.
      alphabet = ['a', 't', 'g']
      1.upto(5) do |length|
        0.upto(3**length - 1) do |mask|
          sequence = length.times.map do |index|
            alphabet[(mask / (3**index)) % 3]
          end.join
          assert_matches_reference(sequence)
        end
      end
    end

    def test_long_homopolymer_is_linearithmic
      length = 200_000
      sequence = Bio::Sequence::NA.new('a' * length)
      elapsed = time_runs(sequence)
      assert_operator(elapsed, :<, 10.0)
      runs = sequence.runs
      assert_equal(1, runs.size)
      assert_equal([0, length, 1], runs.first.to_a)
    end

    def test_long_alternating_is_linearithmic
      length = 200_000
      sequence = Bio::Sequence::NA.new(
        Array.new(length) { |index| index.even? ? 'a' : 't' }.join)
      elapsed = time_runs(sequence)
      assert_operator(elapsed, :<, 10.0)
      runs = sequence.runs
      assert_equal(1, runs.size)
      assert_equal(2, runs.first.period)
      assert_equal([0, length, 2], runs.first.to_a)
    end

    private

    def time_runs(sequence)
      starting = Process.clock_gettime(Process::CLOCK_MONOTONIC)
      sequence.runs
      Process.clock_gettime(Process::CLOCK_MONOTONIC) - starting
    end
  end
end
