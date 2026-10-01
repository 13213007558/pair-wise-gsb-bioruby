# frozen_string_literal: true

# bio/util/rna_fold/partition_function.rb - McCaskill-style partition function
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

require 'bio/util/rna_fold/constants'
require 'bio/util/rna_fold/inside_dp'
require 'bio/util/rna_fold/outside_dp'
require 'bio/util/rna_fold/probabilities'

module Bio
  class RNAFold
    # Log-domain partition function and base-pair edge probabilities for
    # pseudoknot-free RNA secondary structures. The recurrence is a McCaskill
    # style inside/outside dynamic program using only caller supplied additive
    # pairing energies.
    module PartitionFunction
      # Raised when constraints leave no valid secondary structure.
      class InfeasibleConstraintsError < ArgumentError; end

      include InsideDP
      include OutsideDP
      include Probabilities

      # @return [Float] natural log of the partition function Z
      def log_partition_function
        compute! unless @log_partition_function
        @log_partition_function
      end

      # @return [Float] partition function Z
      def partition_function
        Math.exp(log_partition_function)
      end

      # @return [Float] probability that bases +left+ and +right+ are paired
      def pair_probability(left, right)
        compute!
        left = left.to_int
        right = right.to_int
        return 0.0 if left == right
        left, right = [left, right].sort
        check_index(left)
        check_index(right)
        pair_log_probability(left, right)
      end

      # @return [Float] probability that +base+ remains unpaired
      def unpaired_probability(base)
        compute!
        index = base.to_int
        check_index(index)
        return 0.0 if constraints.forced_paired?(index)
        normalized_probability(@log_unpaired[index])
      end

      # @return [Float] probability that +base+ pairs with any other base
      def paired_probability(base)
        1.0 - unpaired_probability(base)
      end

      # @return [Hash{Array(Integer, Integer) => Float}] nonzero pair probs
      def base_pair_probabilities
        compute!
        result = {}
        length.times do |left|
          (left + 1).upto(length - 1) do |right|
            probability = pair_log_probability(left, right)
            result[[left, right].freeze] = probability if probability > 0.0
          end
        end
        result
      end

      private

      def compute!
        @inside = compute_inside
        @log_partition_function = length.zero? ? 0.0 : @inside[0][length - 1]
        unless @log_partition_function > NEG_INFINITY
          raise InfeasibleConstraintsError,
                'constraints do not permit any secondary structure'
        end
        @outside, @enclosing_pair, @enclosing_unpaired =
          compute_outside unless length.zero?
        @log_unpaired = Array.new(length) { |base| compute_log_unpaired(base) }
        self
      end

      def candidate_log_weight(left, right)
        log_weight[left][right]
      end

      def log_sum_exp(terms)
        self.class.log_sum_exp(terms)
      end

      def log_weight
        @log_weight ||= build_log_weights
      end
    end
  end
end
