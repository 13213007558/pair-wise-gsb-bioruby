# frozen_string_literal: true

module Bio
  module DCJ
    class DistanceCalculator
      def initialize(source, target)
        @source = source
        @target = target
      end

      def call
        remaining = Endpoint.all(source.block_ids).to_set
        total = 0
        until remaining.empty?
          start = remaining.first
          total += component_cost(start, remaining)
        end
        total
      end

      private

      attr_reader :source, :target

      def component_cost(start, remaining)
        vertices = component_vertices(start, remaining)
        source_edges = colored_edges(source, vertices)
        target_edges = colored_edges(target, vertices)
        return source_edges - 1 if cycle?(vertices, source_edges, target_edges)

        [source_edges, target_edges].max
      end

      def component_vertices(start, remaining)
        stack = [start]
        remaining.delete(start)
        [start].tap do |vertices|
          visit_stack(stack, vertices, remaining)
        end
      end

      def visit_stack(stack, vertices, remaining)
        until stack.empty?
          point = stack.pop
          collect_neighbors(source, point, remaining, stack, vertices)
          collect_neighbors(target, point, remaining, stack, vertices)
        end
      end

      def collect_neighbors(genome, point, remaining, stack, vertices)
        genome.neighbors(point, remaining).each do |neighbor|
          vertices << neighbor
          stack << neighbor
        end
      end

      def colored_edges(genome, vertices)
        vertex_set = vertices.to_set
        genome.adjacencies.count do |edge|
          edge.all? { |point| vertex_set.include?(point) }
        end
      end

      def cycle?(vertices, source_edges, target_edges)
        source_edges.positive? && source_edges == target_edges &&
          vertices.all? { |point| source.degree(point) == 1 } &&
          vertices.all? { |point| target.degree(point) == 1 }
      end
    end
  end
end
