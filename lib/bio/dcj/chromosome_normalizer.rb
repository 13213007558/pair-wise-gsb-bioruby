# frozen_string_literal: true

module Bio
  module DCJ
    class ChromosomeNormalizer
      def self.linear(blocks)
        orient_minimum_block(blocks)
      end

      def self.circular(blocks)
        oriented = orient_minimum_block(blocks)
        index = oriented.index(oriented.map(&:abs).min)
        oriented.rotate(index)
      end

      def self.orient_minimum_block(blocks)
        minimum = blocks.map(&:abs).min
        reversed = reverse_complement(blocks)
        return blocks if blocks.include?(minimum)
        return reversed if reversed.include?(minimum)

        [blocks, reversed].min
      end

      def self.reverse_complement(blocks)
        blocks.reverse.map(&:-@)
      end
    end
  end
end
