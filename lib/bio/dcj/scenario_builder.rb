# frozen_string_literal: true

module Bio
  module DCJ
    class ScenarioBuilder
      def initialize(source, target)
        @source = source
        @target = target
      end

      def call
        operations = []
        current = source
        until current.adjacencies == target.adjacencies
          operation = next_operation(current)
          operations << operation
          current = current.apply(operation)
        end
        operations
      end

      private

      attr_reader :source, :target

      def next_operation(current)
        current_distance = DCJ.distance(current, target)
        candidate = OperationCandidates.new(current).to_a.find do |operation|
          DCJ.distance(current.apply(operation), target) < current_distance
        end
        candidate || raise('unable to construct a DCJ scenario')
      end
    end
  end
end
