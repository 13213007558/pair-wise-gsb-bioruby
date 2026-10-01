# frozen_string_literal: true

# bio/util/rna_fold/inside_dp.rb - Log-domain inside dynamic programming
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

module Bio
  class RNAFold
    # Computes the inside partition function for one subsequence interval.
    #
    # The grammar is unambiguous and only generates pseudoknot-free,
    # non-crossing structures:
    #   F(i, j) = F(i+1, j)                               (i unpaired)
    #           + sum_{q, pair(i,q), q<=j} w(i,q) F(i+1,q-1) F(q+1,j)
    module InsideDP
      private

      # @return [Array<Array<Float>>] log inside weights; nil means -infinity
      def compute_inside
        table = Array.new(length) { Array.new(length, NEG_INFINITY) }
        fill_inside(table)
        table
      end

      def fill_inside(table)
        (length - 1).downto(0) do |left|
          left.upto(length - 1) do |right|
            table[left][right] = inside_value(left, right, table)
          end
        end
      end

      def inside_value(left, right, table)
        terms = []
        terms << inside_value_unpaired(left, right, table)
        terms.concat(inside_pair_terms(left, right, table))
        log_sum_exp(terms)
      end

      def inside_value_unpaired(left, right, table)
        return NEG_INFINITY if constraints.forced_partner(left)
        inside_at(left + 1, right, table)
      end

      def inside_pair_terms(left, right, table)
        minimum = left + min_hairpin_length + 1
        minimum.upto(right).filter_map do |partner|
          weight = candidate_log_weight(left, partner)
          next unless weight
          weight +
            inside_at(left + 1, partner - 1, table) +
            inside_at(partner + 1, right, table)
        end
      end

      def inside_at(left, right, table)
        return 0.0 if left > right
        table[left][right]
      end
    end
  end
end
