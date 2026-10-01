# frozen_string_literal: true

# bio/codon_optimizer/validator.rb - input validation for codon optimizer

require 'set'

module Bio
  class CodonOptimizer
    # Stateless validation helpers for CodonOptimizer.
    module Validator
      module_function

      def normalize_nucleotides(value, description)
        unless value.respond_to?(:to_str)
          raise ArgumentError, "#{description} must be a String"
        end

        normalized = value.to_str.downcase.gsub(/\s+/, '').tr('u', 't')
        unless normalized.match?(/\A[acgt]*\z/)
          raise ArgumentError, "#{description} contains an invalid nucleotide"
        end

        normalized
      end

      def normalize_patterns(patterns)
        values = patterns.is_a?(String) ? [patterns] : Array(patterns)
        normalized = values.map do |pattern|
          normalize_nucleotides(pattern, 'forbidden substring')
        end.uniq
        normalized.each do |pattern|
          next unless pattern.empty?

          raise ArgumentError, 'forbidden substring cannot be empty'
        end
      end

      def normalize_weights(weights)
        unless weights.respond_to?(:each_pair)
          raise ArgumentError, 'codon_weights must respond to each_pair'
        end

      weights.each_pair.with_object({}) do |(codon, weight), normalized|
        normalized_codon = normalize_nucleotides(codon, 'codon weight')
        unless normalized_codon.length == 3 && weight.is_a?(Numeric)
          raise ArgumentError, "invalid codon weight: #{codon.inspect}"
        end
          normalized[normalized_codon] = weight
        end
      end

      def normalize_positions(positions)
        return Set.new if positions.nil?

      positions = [positions] if positions.is_a?(Integer)
      unless positions.respond_to?(:each)
        raise ArgumentError, 'fixed_positions must be enumerable'
        end

      positions.each_with_object(Set.new) do |position, set|
        Array(position).each do |item|
          validate_position(item)
          set << item
        end
        end
      end

      def normalize_gc_range(range, minimum, maximum)
        if range
          unless range.respond_to?(:min) && range.respond_to?(:max)
            raise ArgumentError, 'gc_count must have min and max methods'
          end
          minimum = range.min
          maximum = range.max
        end

        minimum ||= 0
        validate_gc_value(minimum, 'min_gc_count')
        validate_gc_value(maximum, 'max_gc_count') if maximum
        if maximum && minimum > maximum
          raise ArgumentError, 'min_gc_count cannot exceed max_gc_count'
        end

        [minimum, maximum]
      end

      def validate_strategy(strategy, name)
        unless %i[preserve optimize].include?(strategy)
          raise ArgumentError, "#{name} must be :preserve or :optimize"
        end
        strategy
      end

      def validate_codons(table, codons)
        codons.each_with_index do |codon, index|
          next unless table[codon].nil?

          raise ArgumentError, "codon #{codon.inspect} at position #{index} " \
                               'is not defined by the genetic code table'
        end
      end

      def validate_position(item)
        unless item.is_a?(Integer) && !item.negative?
          raise ArgumentError, "invalid fixed codon position: #{item.inspect}"
        end
        true
      end

      def validate_gc_value(value, name)
        return if value.is_a?(Integer) && !value.negative?

        raise ArgumentError, "invalid #{name}: #{value.inspect}"
      end

    end
  end
end
