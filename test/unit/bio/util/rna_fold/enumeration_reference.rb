# frozen_string_literal: true

# test/unit/bio/util/rna_fold/enumeration_reference.rb - Reference enumerator
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

module Bio
  # Test support that enumerates every pseudoknot-free matching.
  module RNAFoldEnumerationReference
    def enumerate_structures(length, min_hairpin, allowed)
      structures = []
      visit = lambda do |position, pairs|
        if position >= length
          structures << pairs.dup
          return
        end
        occupied = pairs.flatten
        if occupied.include?(position)
          visit.call(position + 1, pairs)
          next
        end
        visit.call(position + 1, pairs)
        add_pairs(position, min_hairpin, allowed, pairs, occupied, visit)
      end
      visit.call(0, [])
      structures
    end

    def add_pairs(position, min_hairpin, allowed, pairs, occupied, visit)
      (position + min_hairpin + 1).upto(allowed.length - 1) do |partner|
        next if occupied.include?(partner) || !allowed[position][partner]
        next if crossing?(pairs, position, partner)
        visit.call(position + 1, pairs + [[position, partner]])
      end
    end

    def crossing?(pairs, position, partner)
      pairs.any? do |left, right|
        (left < position && position < right && right < partner) ||
          (position < left && left < partner && partner < right)
      end
    end

    def allowed_pairs(sequence, min_hairpin, energies)
      allowed = Array.new(sequence.length) { Array.new(sequence.length, false) }
      sequence.length.times do |left|
        (left + 1).upto(sequence.length - 1) do |right|
          next unless right - left - 1 >= min_hairpin
          value = energies.call(left, right, sequence[left], sequence[right])
          allowed[left][right] = !value.nil?
        end
      end
      allowed
    end

    def reference_probabilities(sequence, min_hairpin, energies, folder)
      allowed = allowed_pairs(sequence, min_hairpin, energies)
      structures = enumerate_structures(sequence.length, min_hairpin, allowed)
      accumulate(structures, sequence, energies, folder)
    end

    def accumulate(structures, sequence, energies, folder)
      weights = {}
      unpaired_weight = Array.new(sequence.length, 0.0)
      partition = structures.sum(0.0) do |pairs|
        next 0.0 unless constraints_match?(folder, pairs)
        weight = structure_weight(pairs, sequence, energies)
        collect_pair_weights(pairs, weights, weight)
        collect_unpaired_weight(pairs, unpaired_weight, sequence.length, weight)
        weight
      end
      [partition, weights, unpaired_weight]
    end

    def collect_pair_weights(pairs, weights, weight)
      pairs.each { |pair| weights[pair] = (weights[pair] || 0.0) + weight }
    end

    def collect_unpaired_weight(pairs, unpaired_weight, length, weight)
      paired = pairs.flatten
      length.times do |base|
        unpaired_weight[base] += weight unless paired.include?(base)
      end
    end

    def structure_weight(pairs, sequence, energies)
      pairs.reduce(1.0) do |product, (left, right)|
        energy = energies.call(left, right, sequence[left], sequence[right])
        product * Math.exp(-energy)
      end
    end

    def constraints_match?(folder, pairs)
      constraints = folder.constraints
      return false unless constraints.forced_pairs.all? { |p| pairs.include?(p) }
      return false if pairs.any? { |pair| constraints.forbidden_pairs.include?(pair) }
      constraints.unpaired_positions.none? { |base| pairs.flatten.include?(base) }
    end

    def assert_matches_enumeration(folder, energies, sequence, min_hairpin)
      partition, pair_weights, unpaired_weight =
        reference_probabilities(sequence, min_hairpin, energies, folder)
      assert_in_delta(partition, folder.partition_function, 1e-9)
      assert_base_probabilities(folder, partition, pair_weights,
                                unpaired_weight, sequence.length)
    end

    def assert_base_probabilities(folder, partition, pair_weights,
                                  unpaired_weight, length)
      length.times do |base|
        expected_paired = pair_marginal(folder, base, length, partition, pair_weights)
        expected_unpaired = unpaired_weight[base] / partition
        assert_unpaired_marginal(folder, base, expected_unpaired)
        assert_in_delta(1.0, expected_unpaired + expected_paired, 1e-10,
                        "enumerated normalization #{base}")
        assert_fold_normalization(folder, base)
      end
    end

    def assert_unpaired_marginal(folder, base, expected_unpaired)
      assert_in_delta(expected_unpaired, folder.unpaired_probability(base),
                      1e-10, "unpaired #{base}")
    end

    def assert_fold_normalization(folder, base)
      total = folder.unpaired_probability(base) + folder.paired_probability(base)
      assert_in_delta(1.0, total, 1e-10, "fold normalization #{base}")
    end

    def pair_marginal(folder, base, length, partition, weights)
      expected_paired = 0.0
      length.times do |other|
        next if base == other
        left, right = [base, other].sort
        expected = (weights[[left, right]] || 0.0) / partition
        assert_in_delta(expected, folder.pair_probability(left, right), 1e-10,
                        "pair #{[left, right]}")
        expected_paired += expected
      end
      expected_paired
    end
  end
end
