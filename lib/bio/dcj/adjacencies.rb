# frozen_string_literal: true

module Bio
  module DCJ
    class Adjacencies
      def self.build(chromosomes)
        chromosomes.each_with_object(Set.new) do |chromosome, result|
          chromosome.blocks.each_cons(2) do |first, second|
            result.add(pair_between(first, second))
          end
          if chromosome.circular?
            result.add(pair_between(chromosome.blocks.last, chromosome.blocks.first))
          end
        end
      end

      def self.pair(first, second)
        [first, second].sort
      end

      def self.pair_between(first, second)
        pair(right_end(first), left_end(second))
      end

      def self.right_end(block)
        Endpoint.new(block.abs, block.positive? ? :h : :t)
      end

      def self.left_end(block)
        Endpoint.new(block.abs, block.positive? ? :t : :h)
      end

      def self.validate(edges, block_ids)
        degrees = Endpoint.all(block_ids).to_h { |point| [point, 0] }
        edges.each { |edge| validate_edge(edge, degrees) }
      end

      def self.validate_edge(edge, degrees)
        unless edge.is_a?(Array) && edge.length == 2 && edge.uniq.length == 2
          raise ArgumentError, 'each adjacency must contain two distinct endpoints'
        end
        unless edge.all? { |point| point.is_a?(Endpoint) && degrees.key?(point) }
          raise ArgumentError, 'adjacency contains an unknown endpoint'
        end

        edge.each { |point| count_endpoint(point, degrees) }
      end

      def self.count_endpoint(point, degrees)
        degrees[point] += 1
        return if degrees[point] <= 1

        raise ArgumentError, 'endpoint appears in more than one adjacency'
      end

      private_class_method :count_endpoint
    end
  end
end
