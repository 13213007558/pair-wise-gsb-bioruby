# frozen_string_literal: true

module Bio
  module DCJ
    class OperationCandidates
      def initialize(genome)
        @genome = genome
      end

      def to_a
        candidates = fission_candidates + fusion_candidates +
                     telomere_rejoin_candidates + double_cut_candidates
        candidates.compact
      end

      private

      attr_reader :genome

      def adjacencies
        genome.adjacencies.to_a.sort
      end

      def fission_candidates
        adjacencies.map { |edge| Operation.new([edge], []) }
      end

      def fusion_candidates
        genome.telomeres.sort.combination(2).map do |pair|
          Operation.new([], [pair.sort])
        end
      end

      def telomere_rejoin_candidates
        adjacencies.flat_map do |edge|
          first, second = edge
          genome.telomeres.flat_map do |telomere|
            [
              rejoin_operation(edge, first, telomere),
              rejoin_operation(edge, second, telomere)
            ]
          end
        end
      end

      def rejoin_operation(edge, endpoint, telomere)
        return nil if endpoint == telomere

        Operation.new([edge], [[endpoint, telomere].sort])
      end

      def double_cut_candidates
        adjacencies.combination(2).flat_map do |first, second|
          reconnections(first, second)
        end
      end

      def reconnections(first, second)
        pairings(first, second).map do |pairs|
          reconnection([first, second], pairs)
        end
      end

      def pairings(first, second)
        a, b = first
        c, d = second
        [[[a, c], [b, d]], [[a, d], [b, c]]]
      end

      def reconnection(cuts, pairs)
        return nil if pairs.any? { |first, second| first == second }

        Operation.new(cuts, pairs.map(&:sort))
      end
    end
  end
end
