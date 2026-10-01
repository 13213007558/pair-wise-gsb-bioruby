# frozen_string_literal: true

# test/unit/bio/util/rna_fold/test_partition_function.rb - DP vs enumeration
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

require_relative 'test_enumeration_helper'

module Bio
  class TestRNAFoldPartitionFunction < Test::Unit::TestCase
    include RNAFoldEnumerationReference
    def test_uniform_short_sequence_matches_exhaustive_enumeration
      sequence = 'AAAAAA'
      min_hairpin = 1
      energies = ->(_left, _right, _a, _b) { -1.0 }
      folder = Bio::RNAFold.new(sequence, min_hairpin_length: min_hairpin,
                                pairing_energy: energies)
      assert_matches_enumeration(folder, energies, sequence, min_hairpin)
    end

    def test_random_energies_match_exhaustive_enumeration
      rng = Random.new(1234)
      1.upto(9) do |length|
        10.times { check_random_case(length, rng) }
      end
    end

    def check_random_case(length, rng)
      sequence = 'A' * length
      min_hairpin = rng.rand(0..3)
      table = Array.new(length) { Array.new(length) }
      length.times do |left|
        (left + 1).upto(length - 1) do |right|
          next unless right - left - 1 >= min_hairpin
          table[left][right] = rng.rand < 0.2 ? nil : rng.rand(-3.0..2.0)
        end
      end
      energies = ->(left, right, _a, _b) { table[left] && table[left][right] }
      folder = Bio::RNAFold.new(sequence, min_hairpin_length: min_hairpin,
                                pairing_energy: energies)
      assert_matches_enumeration(folder, energies, sequence, min_hairpin)
    end

    def test_constraints_match_exhaustive_enumeration
      sequence = 'AAAAAAA'
      min_hairpin = 1
      energies = ->(_left, _right, _a, _b) { -0.7 }
      constraint_cases.each do |options|
        folder = Bio::RNAFold.new(sequence, min_hairpin_length: min_hairpin,
                                  pairing_energy: energies, **options)
        assert_matches_enumeration(folder, energies, sequence, min_hairpin)
      end
    end

    def constraint_cases
      [
        { forced_pairs: [[0, 6]] },
        { forced_pairs: [[0, 6], [2, 4]] },
        { forbidden_pairs: [[0, 6], [1, 6]] },
        { unpaired: [0, 3, 6] },
        { forced_pairs: [[0, 6]], unpaired: [3], forbidden_pairs: [[1, 5]] }
      ]
    end

    def test_minimum_hairpin_length_is_enforced
      energies = ->(_left, _right, _a, _b) { -100.0 }
      folder = Bio::RNAFold.new('AAAA', min_hairpin_length: 3,
                                pairing_energy: energies)
      assert_in_delta(1.0, folder.partition_function, 1e-12)
      assert_in_delta(1.0, folder.unpaired_probability(0), 1e-12)
    end

    def test_beta_scales_pair_weights
      sequence = 'AAA'
      energies = ->(_left, _right, _a, _b) { -1.0 }
      weak = new_folder(sequence, energies, 0.5)
      strong = Bio::RNAFold.new(sequence, min_hairpin_length: 0,
                                pairing_energy: energies, beta: 2.0)
      assert weak.pair_probability(0, 2) < strong.pair_probability(0, 2)
      assert_in_delta(1.0 + 3.0 * Math.exp(0.5), weak.partition_function,
                      1e-12)
    end

    def new_folder(sequence, energies, beta)
      Bio::RNAFold.new(sequence, min_hairpin_length: 0,
                       pairing_energy: energies, beta: beta)
    end

    def test_hash_and_matrix_energy_sources_agree
      hash_energy = { [0, 3] => -1.5, [1, 2] => 0.5 }
      from_hash = Bio::RNAFold.new('AAAA', min_hairpin_length: 0,
                                   pairing_energy: hash_energy)
      from_matrix = Bio::RNAFold.new('AAAA', min_hairpin_length: 0,
                                    pairing_energy: energy_matrix)
      assert_in_delta(from_hash.partition_function, from_matrix.partition_function, 1e-12)
      assert_in_delta(from_hash.pair_probability(0, 3), from_matrix.pair_probability(3, 0), 1e-12)
    end

    def energy_matrix
      matrix = Array.new(4) { Array.new(4) }
      matrix[0][3] = -1.5
      matrix[1][2] = 0.5
      matrix
    end
  end
end
