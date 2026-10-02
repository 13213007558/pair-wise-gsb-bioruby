# frozen_string_literal: true

module Bio
  module DCJ
    class Genome
      attr_reader :chromosomes, :adjacencies, :block_ids, :telomeres

      def initialize(raw_chromosomes)
        parsed = GenomeParser.new(raw_chromosomes).parse
        @block_ids = GenomeParser.block_ids(parsed)
        @adjacencies = Adjacencies.build(parsed).freeze
        initialize_derived_state
      end

      def apply(operation)
        cuts, joins = OperationValidator.new(
          adjacencies, telomeres, operation
        ).call
        next_adjacencies = adjacencies.dup
        cuts.each { |edge| next_adjacencies.delete(edge) }
        joins.each { |edge| next_adjacencies.add(edge) }
        self.class.from_adjacencies(next_adjacencies, block_ids)
      end

      def self.from_adjacencies(adjacencies, block_ids)
        allocate.tap do |genome|
          genome.send(:initialize_from_adjacencies, adjacencies.to_set, block_ids)
        end
      end

      def endpoint_adjacencies
        adjacencies.to_a
      end

      def neighbors(point, remaining)
        endpoint_adjacencies.filter_map do |first, second|
          neighbor = first == point ? second : second == point ? first : nil
          neighbor if neighbor && remaining.delete?(neighbor)
        end
      end

      def degree(point)
        endpoint_adjacencies.count { |edge| edge.include?(point) }
      end

      def ==(other)
        other.is_a?(Genome) && adjacencies == other.adjacencies &&
          block_ids == other.block_ids
      end
      alias eql? ==

      def hash
        [block_ids, adjacencies].hash
      end

      private

      def initialize_from_adjacencies(edges, ids)
        @block_ids = ids.uniq.sort.freeze
        Adjacencies.validate(edges, block_ids)
        @adjacencies = edges.freeze
        initialize_derived_state
      end

      def initialize_derived_state
        endpoints = Endpoint.all(block_ids)
        @telomeres = endpoints.select { |point| degree(point).zero? }.sort.freeze
        @chromosomes = ChromosomeBuilder.new(adjacencies, block_ids).build.freeze
        @block_ids = block_ids.freeze
      end
    end
  end
end
