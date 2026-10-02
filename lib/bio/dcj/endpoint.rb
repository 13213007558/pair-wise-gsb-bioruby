# frozen_string_literal: true

module Bio
  module DCJ
    Endpoint = Struct.new(:block, :side) do
      include Comparable

      def <=>(other)
        return nil unless other.is_a?(Endpoint)

        [block, side] <=> [other.block, other.side]
      end

      def opposite
        Endpoint.new(block, side == :h ? :t : :h)
      end

      def self.all(block_ids)
        block_ids.flat_map { |id| [new(id, :h), new(id, :t)] }
      end
    end
  end
end
