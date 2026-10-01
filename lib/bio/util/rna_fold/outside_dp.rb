# frozen_string_literal: true

# bio/util/rna_fold/outside_dp.rb - Log-domain outside dynamic programming
#
# Copyright:: Copyright (c) 2026 BioRuby Project
# License:: The Ruby License

module Bio
  class RNAFold
    # Computes log outside weights +O(i, j)+ for every possible pair.
    #
    # +O(i, j)+ is the total Boltzmann weight of the part outside pair
    # +(i, j)+, including every parent pair and every sibling of the pair.
    # Pair edge probability is therefore
    # +exp(w(i,j) + inside(i+1,j-1) + O(i,j) - log_Z)+.
    #
    # Rows are processed by increasing left endpoint. For each parent left
    # +h+ two accumulators collect parents while their right endpoint +l+
    # grows, making both edge and unpaired probabilities O(n**3).
    module OutsideDP
      private

      def compute_outside
        outside = Array.new(length) { Array.new(length, NEG_INFINITY) }
        enclosing_pair = Array.new(length) { Array.new(length, NEG_INFINITY) }
        enclosing_unpaired = Array.new(length) { Array.new(length, NEG_INFINITY) }
        fill_outside(outside, enclosing_pair, enclosing_unpaired)
        [outside, enclosing_pair, enclosing_unpaired]
      end

      def fill_outside(outside, enclosing_pair, enclosing_unpaired)
        length.times do |left|
          fill_outside_row(left, outside, enclosing_pair)
          update_enclosing(left, outside, enclosing_pair, enclosing_unpaired)
        end
      end

      def fill_outside_row(left, outside, enclosing_pair)
        (left + 1).upto(length - 1) do |right|
          next unless candidate_log_weight(left, right)
          outside[left][right] = log_sum_exp([
            top_outside(left, right),
            enclosing_pair_terms(left, right, enclosing_pair)
          ])
        end
      end

      def enclosing_pair_terms(child_left, child_right, enclosing)
        terms = (0...child_left).filter_map do |parent_left|
          value = enclosing[parent_left][child_right]
          next if value == NEG_INFINITY
          inside_at_value(parent_left + 1, child_left - 1) + value
        end
        log_sum_exp(terms)
      end

      def top_outside(left, right)
        inside_prefix(left) + inside_suffix(right)
      end

      def inside_prefix(left)
        left.zero? ? 0.0 : @inside[0][left - 1]
      end

      def inside_suffix(right)
        right == length - 1 ? 0.0 : @inside[right + 1][length - 1]
      end

      def update_enclosing(left, outside, enclosing_pair, enclosing_unpaired)
        (left + 1).upto(length - 1) do |right|
          next unless candidate_log_weight(left, right)
          base = outside[left][right] + candidate_log_weight(left, right)
          add_pair_contributions(left, right, base, enclosing_pair)
          add_unpaired_contributions(left, right, base, enclosing_unpaired)
        end
      end

      def add_pair_contributions(parent_left, parent_right, base, enclosing)
        (parent_left + 1).upto(parent_right - 1) do |child_right|
          value = base + inside_at_value(child_right + 1, parent_right - 1)
          add_to_enclosing(enclosing, parent_left, child_right, value)
        end
      end

      def add_unpaired_contributions(parent_left, parent_right, base,
                                     enclosing)
        (parent_left + 1).upto(parent_right - 1) do |child_left|
          next if constraints.forced_partner(child_left)
          value = base + inside_at_value(child_left + 1, parent_right - 1)
          add_to_enclosing(enclosing, parent_left, child_left, value)
        end
      end

      def add_to_enclosing(enclosing, parent_left, child, value)
        enclosing[parent_left][child] = log_sum_exp([
          enclosing[parent_left][child], value
        ])
      end

      def inside_at_value(left, right)
        return 0.0 if left > right
        @inside[left][right]
      end
    end
  end
end
