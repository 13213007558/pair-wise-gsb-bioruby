# frozen_string_literal: true

module Bio
  module DCJ
    class ChromosomeBuilder
      def initialize(adjacencies, block_ids)
        @adjacencies = adjacencies
        @block_ids = block_ids
        @traversal = ChromosomeTraversal.new(adjacencies)
      end

      def build
        remaining = Endpoint.all(block_ids).to_set
        build_linear(remaining) + build_circular(remaining)
      end

      private

      attr_reader :adjacencies, :block_ids, :traversal

      def build_linear(remaining)
        telomeres.filter_map do |point|
          next unless remaining.delete?(point)

          blocks = traversal.linear(point, remaining)
          Chromosome.new(ChromosomeNormalizer.linear(blocks), false)
        end
      end

      def build_circular(remaining)
        Endpoint.all(block_ids).filter_map do |point|
          next unless remaining.delete?(point)

          blocks = traversal.circular(point, remaining)
          Chromosome.new(ChromosomeNormalizer.circular(blocks), true)
        end
      end

      def telomeres
        Endpoint.all(block_ids).select do |point|
          adjacencies.none? { |edge| edge.include?(point) }
        end
      end
    end
  end
end
