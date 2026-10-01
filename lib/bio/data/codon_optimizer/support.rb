# frozen_string_literal: true

module Bio
  class CodonOptimizer
    # Validation helpers for caller-provided optimizer arguments.
    module Arguments

      module_function

      def weights(values)
        return {} if values.nil?

        unless values.is_a?(Hash)
          raise ArgumentError, 'codon_weights must be a Hash'
        end

        values.each_with_object({}) do |(codon, weight), normalized|
          valid_weight(codon, weight)
          normalized[codon.downcase.tr('u', 't')] = weight
        end
      end

      def patterns(values)
        assert_array(values, 'forbidden_substrings')
        values.map { |value| pattern(value) }.uniq
      end

      def positions(values)
        assert_array(values, 'fixed_positions')
        expanded = values.flat_map { |value| expand_position(value) }
        if expanded.any?(&:negative?)
          raise ArgumentError, 'fixed positions must not be negative'
        end

        expanded.uniq.sort
      end

      def gc_range(range)
        return nil if range.nil?

        unless integer_range?(range)
          raise ArgumentError, 'gc_range must be an Integer Range'
        end
        if range.begin.negative? || range.begin > range.end
          raise ArgumentError, "invalid gc_range: #{range.inspect}"
        end

        range
      end

      def valid_weight(codon, weight)
        valid_codon(codon)
        return if weight.is_a?(Numeric)

        raise ArgumentError, "weight for #{codon} is not numeric"
      end

      def valid_codon(codon)
        return if codon.is_a?(String) &&
                  codon.length == 3 &&
                  codon.downcase.match?(/\A[atgc]{3}\z/)

        raise ArgumentError, "invalid codon in codon_weights: #{codon.inspect}"
      end

      def assert_array(value, name)
        return if value.is_a?(Array)

        raise ArgumentError, "#{name} must be an Array"
      end

      def pattern(value)
        unless value.is_a?(String) && !value.empty?
          raise ArgumentError, 'forbidden substring must be a nonempty String'
        end

        normalized = value.downcase.tr('u', 't')
        return normalized if normalized.match?(/\A[atgc]+\z/)

        raise ArgumentError, "invalid forbidden substring: #{value}"
      end

      def expand_position(value)
        if value.is_a?(Integer)
          [value]
        elsif integer_range?(value)
          value.to_a
        else
          raise ArgumentError, "invalid fixed position: #{value.inspect}"
        end
      end

      def integer_range?(value)
        value.is_a?(Range) &&
          value.begin.is_a?(Integer) &&
          value.end.is_a?(Integer)
      end

    end

    # Shared conflict, result, and nucleotide-count helpers.
    module Support

      module_function

      def conflict(type, message, **details)
        details.merge(type: type, message: message)
      end

      def gc_count(sequence)
        sequence.count('g') + sequence.count('c')
      end

      def feasible(result_class, result, rna)
        sequence = result.fetch(:codons).join
        sequence.tr!('t', 'u') if rna
        result_class.new(
          sequence: sequence,
          score: result.fetch(:score),
          gc_count: result.fetch(:gc_count),
          conflicts: []
        )
      end

      def infeasible(result_class, conflicts)
        result_class.new(sequence: nil, score: nil, gc_count: nil,
                         conflicts: conflicts)
      end

    end

  end

end
