# frozen_string_literal: true

# bio/util/rna_fold/energy_adapter.rb - Pairing energy source for RNA folding
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

module Bio
  class RNAFold
    # Adapts caller supplied pairing energies to a uniform +call(i, j)+ API.
    #
    # Supported energy sources are:
    # * +nil+ - every pair is forbidden
    # * a callable object called with +(left, right, left_base, right_base)+
    # * a +Hash+ keyed by +[left, right]+ (keys are orientation independent)
    # * an indexable matrix, such as an +Array+ of arrays
    #
    # Returning +nil+ means that the pair is not allowed. Finite energies are
    # additive Boltzmann energies; negative values are stabilizing. NaN and
    # infinite values are rejected eagerly rather than poisoning probabilities.
    class EnergyAdapter
      # @param source [nil, #call, #[], Hash] pairing energy source
      # @param sequence [String] sequence for callable based energies
      def initialize(source, sequence)
        @sequence = sequence
        @source = source
        @lookup = build_lookup
      end

      # @return [Numeric, nil] finite additive pairing energy, or nil
      def call(left, right)
        value = @lookup.call(left, right)
        return nil if value.nil?
        unless value.is_a?(Numeric)
          raise TypeError, "pairing energy must be numeric, got #{value.class}"
        end
        validate_value(value)
        value
      end

      private

      def build_lookup
        case @source
        when nil then method(:nil_energy)
        when Hash then method(:hash_energy)
        else
          if @source.respond_to?(:call)
            method(:callable_energy)
          else
            method(:matrix_energy)
          end
        end
      end

      def nil_energy(_left, _right)
        nil
      end

      def callable_energy(left, right)
        @source.call(left, right, @sequence[left], @sequence[right])
      end

      def hash_energy(left, right)
        return @source[[left, right]] if @source.key?([left, right])
        @source.key?([right, left]) ? @source[[right, left]] : nil
      end

      def matrix_energy(left, right)
        @source[left][right]
      rescue NoMethodError, TypeError
        raise TypeError, "unsupported pairing energy source: #{@source.class}"
      end

      def validate_value(value)
        return if value.finite?
        raise ArgumentError, "pairing energy must be finite, got #{value}"
      end
    end
  end
end
