# frozen_string_literal: true

# Parameter validation helpers used by Bio::PCR.
module Bio
  class PCR
    # Validation for PCR simulation parameters.
    module Options
      module_function

      def validate(max_mismatches, max_total_mismatches,
                   three_prime_mismatches, three_prime_region,
                   min_product_length, max_product_length)
        non_negative_integer(:max_mismatches, max_mismatches)
        unless max_total_mismatches.nil?
          non_negative_integer(:max_total_mismatches,
                               max_total_mismatches)
        end
        non_negative_integer(:three_prime_mismatches,
                             three_prime_mismatches)
        non_negative_integer(:three_prime_region, three_prime_region)
        positive_integer(:min_product_length, min_product_length)
        unless max_product_length.nil?
          positive_integer(:max_product_length, max_product_length)
        end
        return if three_prime_mismatches <= max_mismatches

        raise ArgumentError,
              'three_prime_mismatches cannot be greater than max_mismatches'
      end

      def validate_sequences(template, forward, reverse, circular,
                             min_product_length, max_product)
        if template.empty? || forward.empty? || reverse.empty?
          raise ArgumentError, 'template and primers must not be empty'
        end

        if circular && (forward.length > template.length ||
                        reverse.length > template.length)
          raise ArgumentError,
                'primers cannot be longer than a circular template'
        end
        return if max_product >= min_product_length

        raise ArgumentError,
              'min_product_length cannot be greater than max_product_length'
      end

      def non_negative_integer(name, value)
        return if value.is_a?(Integer) && value >= 0

        raise ArgumentError, "#{name} must be a non-negative integer"
      end

      def positive_integer(name, value)
        return if value.is_a?(Integer) && value > 0

        raise ArgumentError, "#{name} must be a positive integer"
      end
    end
  end
end
