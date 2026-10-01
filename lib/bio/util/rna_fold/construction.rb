# frozen_string_literal: true

# bio/util/rna_fold/construction.rb - Construction helpers for RNAFold
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

module Bio
  class RNAFold
    # Private initialization helpers shared by the RNAFold constructor.
    module Construction
      private

      def convert_sequence(sequence)
        value = sequence.respond_to?(:to_str) ? sequence.to_str : sequence.to_s
        value.dup.freeze
      rescue NoMethodError
        raise TypeError,
              "sequence must be convertible to a String, got #{sequence.class}"
      end

      def validate_energies!
        build_log_weights
      end

      def prepare_energies(pairing_energy)
        @energy = EnergyAdapter.new(pairing_energy, @sequence)
        validate_energies!
        validate_forced_pairs!
      end

      def build_log_weights
        table = Array.new(length) { Array.new(length) }
        each_candidate_pair do |left, right, energy|
          table[left][right] = -@beta * energy
        end
        table
      end

      def each_candidate_pair
        length.times do |left|
          (left + min_hairpin_length + 1).upto(length - 1) do |right|
            next unless pair_allowed?(left, right)
            energy = energy_at(left, right)
            next if energy.nil?
            yield(left, right, energy)
          end
        end
      end

      def validate_forced_pairs!
        constraints.forced_pairs.each do |left, right|
          validate_forced_pair(left, right)
        end
      end

      def validate_forced_pair(left, right)
        unless right - left - 1 >= min_hairpin_length
          raise Constraints::ConflictError,
                "forced pair (#{left}, #{right}) violates minimum hairpin length"
        end
        return if pair_allowed?(left, right) && !energy_at(left, right).nil?
        raise Constraints::ConflictError,
              "forced pair (#{left}, #{right}) is not a permitted pair"
      end

      def energy_at(left, right)
        @energy.call(left, right)
      end

      def pair_allowed?(left, right)
        constraints.pair_allowed?(left, right)
      end
    end
  end
end
