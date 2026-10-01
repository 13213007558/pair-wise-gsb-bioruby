# frozen_string_literal: true

# bio/util/rna_fold/probabilities.rb - Edge probability extraction
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

module Bio
  class RNAFold
    # Converts finished inside/outside tables into base-pair and unpaired
    # marginal probabilities.
    module Probabilities
      private

      def pair_log_probability(left, right)
        weight = candidate_log_weight(left, right)
        return 0.0 unless weight
        log_weight = weight + inside_at_value(left + 1, right - 1) +
                     @outside[left][right] - @log_partition_function
        normalized_probability(log_weight)
      end

      def compute_log_unpaired(base)
        terms = [unpaired_top_term(base)]
        terms.concat(unpaired_parent_terms(base))
        log_sum_exp(terms) - @log_partition_function
      end

      # Base is the leftmost unpaired item of a top-level decomposition.
      def unpaired_top_term(base)
        return NEG_INFINITY if constraints.forced_paired?(base)
        prefix = base.zero? ? 0.0 : @inside[0][base - 1]
        suffix = base == length - 1 ? 0.0 : @inside[base + 1][length - 1]
        prefix + suffix
      end

      # Base is unpaired inside an enclosing pair (parent_left, right).
      def unpaired_parent_terms(base)
        (0...base).filter_map do |parent_left|
          value = @enclosing_unpaired[parent_left][base]
          next if value == NEG_INFINITY
          inside_at_value(parent_left + 1, base - 1) + value
        end
      end

      def normalized_probability(log_probability)
        return 0.0 if log_probability == NEG_INFINITY
        probability = Math.exp(log_probability)
        if probability > 1.0 && probability - 1.0 < 1e-12
          return 1.0
        end
        [[probability, 0.0].max, 1.0].min
      end

      def inside_at_value(left, right)
        return 0.0 if left > right
        @inside[left][right]
      end
    end
  end
end
