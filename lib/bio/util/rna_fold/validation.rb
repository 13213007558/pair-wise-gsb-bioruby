# frozen_string_literal: true

# bio/util/rna_fold/validation.rb - Numeric validation helpers
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

module Bio
  class RNAFold
    # Private index and parameter validation helpers.
    module Validation
      private

      def check_index(index)
        return if valid_index?(index)
        raise IndexError,
              "base index #{index} outside sequence of length #{length}"
      end

      def valid_index?(index)
        index >= 0 && index < length
      end

      def validate_nonnegative_integer(value, name)
        integer = value.to_int
        return integer if integer >= 0
        raise ArgumentError, "#{name} must be non-negative"
      rescue NoMethodError
        raise TypeError, "#{name} must be an Integer"
      end

      def validate_positive_number(value, name)
        unless value.is_a?(Numeric) && value.finite? && value > 0
          raise ArgumentError, "#{name} must be a positive finite number"
        end
        value
      end
    end
  end
end
