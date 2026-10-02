# frozen_string_literal: true

module Bio
  module DCJ
    class OperationValidator
      def initialize(adjacencies, telomeres, operation)
        @adjacencies = adjacencies
        @telomeres = telomeres
        @operation = operation
      end

      def call
        cuts = normalize_edges(operation.cuts, 'cuts')
        joins = normalize_edges(operation.joins, 'joins')
        validate_current_edges(cuts)
        validate_new_edges(cuts, joins)
        validate_shape(cuts, joins)
        [cuts, joins]
      end

      private

      attr_reader :adjacencies, :telomeres, :operation

      def validate_current_edges(cuts)
        return if cuts.all? { |edge| adjacencies.include?(edge) }

        raise ArgumentError, 'operation cuts must be current adjacencies'
      end

      def validate_new_edges(cuts, joins)
        invalid = joins.any? { |edge| adjacencies.include?(edge) && !cuts.include?(edge) }
        raise ArgumentError, 'operation joins an existing adjacency' if invalid
      end

      def validate_shape(cuts, joins)
        joined = joins.flatten
        validate_join_endpoints(joined)
        validate_freed_endpoints(cuts, joined)
        return if valid_shape?(cuts.length, joins.length)

        raise ArgumentError, 'invalid DCJ operation shape'
      end

      def validate_join_endpoints(joined)
        return if joined.each_slice(2).all? { |edge| edge.uniq.length == 2 }

        raise ArgumentError, 'join endpoints must be distinct'
      end

      def validate_freed_endpoints(cuts, joined)
        freed = cuts.flatten + telomeres
        return if joined.all? { |point| freed.include?(point) }

        raise ArgumentError, 'joined endpoints must be cut or telomeric'
      end

      def valid_shape?(cut_count, join_count)
        (cut_count == 0 && join_count == 1) ||
          (cut_count == 1 && [0, 1].include?(join_count)) ||
          (cut_count == 2 && join_count == 2)
      end

      def normalize_edges(edges, name)
        raise ArgumentError, "#{name} must be an Array" unless edges.is_a?(Array)

        edges.map { |edge| normalize_edge(edge, name) }
      end

      def normalize_edge(edge, name)
        valid = edge.is_a?(Array) && edge.length == 2 &&
                edge.all? { |point| point.is_a?(Endpoint) }
        raise ArgumentError, "each #{name[0..-2]} must contain two endpoints" unless valid

        edge.sort
      end
    end
  end
end
