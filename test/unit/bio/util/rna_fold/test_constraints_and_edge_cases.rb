# frozen_string_literal: true

# test/unit/bio/util/rna_fold/test_constraints_and_edge_cases.rb
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 4,
                            'bioruby_test_helper.rb')).cleanpath.to_s

require 'test/unit'
require 'bio/util/rna_fold'

module Bio
  class TestRNAFoldConstraints < Test::Unit::TestCase
    def test_conflicting_constraints_are_rejected
      energies = ->(_left, _right, _a, _b) { 0.0 }
      assert_crossing_forced_pairs_rejected(energies)
      assert_forced_and_forbidden_rejected(energies)
      assert_forced_and_unpaired_rejected(energies)
      assert_forced_pair_without_energy_rejected
    end

    def assert_crossing_forced_pairs_rejected(energies)
      assert_raise(Bio::RNAFold::Constraints::ConflictError) do
        Bio::RNAFold.new('AAAAAA', min_hairpin_length: 1,
                         pairing_energy: energies,
                         forced_pairs: [[0, 3], [1, 4]])
      end
    end

    def assert_forced_and_forbidden_rejected(energies)
      assert_raise(Bio::RNAFold::Constraints::ConflictError) do
        Bio::RNAFold.new('AAAAAA', min_hairpin_length: 1,
                         pairing_energy: energies,
                         forced_pairs: [[0, 3]], forbidden_pairs: [[3, 0]])
      end
    end

    def assert_forced_and_unpaired_rejected(energies)
      assert_raise(Bio::RNAFold::Constraints::ConflictError) do
        Bio::RNAFold.new('AAAAAA', min_hairpin_length: 1,
                         pairing_energy: energies,
                         forced_pairs: [[0, 3]], unpaired: [0])
      end
    end

    def assert_forced_pair_without_energy_rejected
      assert_raise(Bio::RNAFold::Constraints::ConflictError) do
        Bio::RNAFold.new('AAAAAA', pairing_energy: nil,
                         forced_pairs: [[0, 3]])
      end
    end

    def test_forced_pair_probabilities_are_one
      energies = ->(_left, _right, _a, _b) { -1.0 }
      folder = Bio::RNAFold.new('AAAAAAA', min_hairpin_length: 1,
                                pairing_energy: energies,
                                forced_pairs: [[0, 6], [2, 4]])
      assert_in_delta(1.0, folder.pair_probability(0, 6), 1e-12)
      assert_in_delta(1.0, folder.pair_probability(6, 0), 1e-12)
      assert_in_delta(1.0, folder.pair_probability(2, 4), 1e-12)
      [0, 2, 4, 6].each { |base| assert_forced_base_marginal(folder, base) }
    end

    def assert_forced_base_marginal(folder, base)
      assert_in_delta(0.0, folder.unpaired_probability(base), 1e-12)
      assert_in_delta(1.0, folder.paired_probability(base), 1e-12)
    end

    def test_hairpin_too_small_forced_pair_is_rejected
      energies = ->(_left, _right, _a, _b) { -100.0 }
      assert_raise(Bio::RNAFold::Constraints::ConflictError) do
        Bio::RNAFold.new('AAAA', min_hairpin_length: 3,
                         pairing_energy: energies,
                         forced_pairs: [[0, 3]])
      end
    end
  end

  class TestRNAFoldEdgeCases < Test::Unit::TestCase
    def test_infeasible_structure_raises_instead_of_returning_zero
      folder = Bio::RNAFold.new('AAA', min_hairpin_length: 0,
                                pairing_energy: ->(*) { -1.0 },
                                forced_pairs: [[0, 2]])
      folder.send(:log_weight)[0][2] = nil
      assert_raise(
        Bio::RNAFold::PartitionFunction::InfeasibleConstraintsError
      ) { folder.partition_function }
    end

    def test_numerically_stable_with_large_energies
      energies = ->(left, right, _a, _b) do
        left + right == 13 ? -300.0 : 300.0
      end
      folder = Bio::RNAFold.new('A' * 14, min_hairpin_length: 3,
                                pairing_energy: energies)
      assert_in_delta(1.0, folder.pair_probability(0, 13), 1e-12)
      14.times do |base|
        total = folder.unpaired_probability(base) +
                folder.paired_probability(base)
        assert_in_delta(1.0, total, 1e-10)
      end
      assert_equal(true, folder.log_partition_function.finite?)
    end

    def test_empty_sequence_has_one_empty_structure
      folder = Bio::RNAFold.new('', pairing_energy: nil)
      assert_in_delta(0.0, folder.log_partition_function, 1e-15)
      assert_in_delta(1.0, folder.partition_function, 1e-15)
      assert_equal({}, folder.base_pair_probabilities)
    end

    def test_invalid_energy_values_are_rejected
      assert_raise(ArgumentError) do
        Bio::RNAFold.new('AAAA', min_hairpin_length: 0,
                         pairing_energy: ->(*) { Float::NAN })
      end
      assert_raise(ArgumentError) do
        Bio::RNAFold.new('AAAA', min_hairpin_length: 0,
                         pairing_energy: ->(*) { -Float::INFINITY })
      end
    end

    def test_invalid_indices_are_rejected
      folder = Bio::RNAFold.new('AAA', min_hairpin_length: 0,
                                pairing_energy: ->(*) { 0.0 })
      assert_raise(IndexError) { folder.pair_probability(0, 3) }
      assert_raise(IndexError) { folder.pair_probability(-1, 1) }
      assert_raise(IndexError) { folder.unpaired_probability(3) }
    end
  end
end
