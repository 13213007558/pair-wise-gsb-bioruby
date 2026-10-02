# frozen_string_literal: true

module Bio
  module DCJ
    class ChromosomeTraversal
      def initialize(adjacencies)
        @adjacencies = adjacencies.to_a.freeze
      end

      def linear(point, remaining)
        blocks = []
        while point
          mark_block(point, remaining)
          blocks << signed_block(point, point.side)
          point, _side = available_neighbor(point.opposite, remaining)
        end
        blocks
      end

      def circular(point, remaining)
        blocks = []
        first = point
        while point
          point = advance_circle(point, first, remaining, blocks)
        end
        blocks
      end

      private

      attr_reader :adjacencies

      def advance_circle(point, first, remaining, blocks)
        mark_block(point, remaining)
        blocks << signed_block(point, point.side)
        neighbor, _side = ordinary_neighbor(point.opposite)
        return nil if neighbor == first

        neighbor if neighbor && remaining.include?(neighbor)
      end

      def mark_block(point, remaining)
        remaining.delete(point)
        remaining.delete(point.opposite)
      end

      def available_neighbor(point, remaining)
        neighbor, side = ordinary_neighbor(point)
        [neighbor, side] if neighbor && remaining.include?(neighbor)
      end

      def ordinary_neighbor(point)
        edge = adjacencies.find { |candidate| candidate.include?(point) }
        return nil unless edge

        first, second = edge
        neighbor = first == point ? second : first
        [neighbor, neighbor.side]
      end

      def signed_block(point, entered_side)
        entered_side == :h ? point.block : -point.block
      end
    end
  end
end
