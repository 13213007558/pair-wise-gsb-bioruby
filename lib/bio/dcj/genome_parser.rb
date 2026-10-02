# frozen_string_literal: true

module Bio
  module DCJ
    class GenomeParser
      def initialize(raw_chromosomes)
        @raw_chromosomes = raw_chromosomes
      end

      def parse
        unless raw_chromosomes.is_a?(Array)
          raise ArgumentError, 'chromosomes must be an Array'
        end

        raw_chromosomes.map { |item| coerce_chromosome(item) }
      end

      def self.block_ids(chromosomes)
        ids = chromosomes.flat_map { |chromosome| chromosome.blocks.map(&:abs) }
        duplicates = ids.group_by(&:itself)
                       .select { |_, group| group.size > 1 }.keys
        raise ArgumentError, "duplicated blocks: #{duplicates.join(', ')}" unless duplicates.empty?

        ids.sort
      end

      private

      attr_reader :raw_chromosomes

      def coerce_chromosome(item)
        if item.is_a?(Chromosome)
          return Chromosome.new(parse_blocks(item.blocks), circular_value(item.circular))
        end
        return parse_array(item) if item.is_a?(Array)
        return parse_hash(item) if item.is_a?(Hash)

        raise ArgumentError, 'each chromosome must be an Array or a Chromosome'
      end

      def parse_array(item)
        Chromosome.new(parse_blocks(item), false)
      end

      def parse_hash(item)
        keys_valid = item.keys == [:blocks] || item.keys.sort == [:blocks, :circular]
        raise ArgumentError, 'invalid chromosome representation' unless keys_valid

        Chromosome.new(parse_blocks(item[:blocks]),
                       circular_value(item.fetch(:circular, false)))
      end

      def circular_value(value)
        return value if value == true || value == false

        raise ArgumentError, 'chromosome circular flag must be true or false'
      end

      def parse_blocks(blocks)
        unless blocks.is_a?(Array) && !blocks.empty?
          raise ArgumentError, 'chromosome blocks must be a non-empty Array'
        end

        blocks.map { |block| parse_block(block) }
      end

      def parse_block(block)
        return block if block.is_a?(Integer) && block.nonzero?

        raise ArgumentError, "invalid block orientation: #{block.inspect}"
      end
    end
  end
end
