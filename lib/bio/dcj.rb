#
# = bio/dcj.rb - double-cut-and-join genome rearrangements
#
# Copyright:: Copyright (C) 2026 BioRuby Project
# License:: The Ruby License
#
# frozen_string_literal: true

require 'set'

module Bio
  # Double-cut-and-join distance and sorting operations for signed gene
  # block genomes containing linear and circular chromosomes.
  class DCJ
    TAIL = :tail
    HEAD = :head

    # Signed occurrence of a gene block in a chromosome.
    SignedBlock = Struct.new(:id, :orientation, keyword_init: true)

    # Linear or circular chromosome represented by signed blocks.
    ChromosomeSpec = Struct.new(:blocks, :circular, keyword_init: true)

    # Genome represented as a set of block-end adjacencies.
    class Genome
      include Enumerable

      attr_reader :chromosomes, :blocks

      def self.[](*chromosomes)
        new(chromosomes)
      end

      def self.from_adjacencies(adjacencies)
        unless adjacencies.is_a?(Array)
          raise ArgumentError, 'adjacencies must be an array'
        end

        edges = adjacencies.map { |edge| normalize_edge(edge) }
        endpoints = edges.flatten(1)
        block_endpoints = endpoints
        blocks = block_endpoints.map(&:first).uniq
        grouped = block_endpoints.group_by(&:first)
        expected_ends = [HEAD, TAIL].sort_by(&:to_s)
        unless grouped.values.all? { |values| values.map(&:last).sort_by(&:to_s) == expected_ends }
          raise ArgumentError, 'each block must have a head and a tail'
        end

        adjacency = {}
        edges.each do |edge|
          edge.each do |endpoint|
            adjacency[endpoint] = edge
          end
        end
        unless adjacency.size == endpoints.size
          raise ArgumentError, 'each block end can occur at most once'
        end

        chromosomes = []
        visited = Set.new

        endpoints.select { |endpoint| adjacency[endpoint].size == 1 }.sort_by(&:inspect).each do |endpoint|
          next if visited.include?(endpoint)

          blocks_in_chromosome = traverse_linear(adjacency, endpoint, visited)
          chromosomes << ChromosomeSpec.new(blocks: blocks_in_chromosome,
                                            circular: false)
        end

        blocks.sort_by(&:inspect).each do |block|
          start = [block, TAIL]
          next if visited.include?(start)

          blocks_in_chromosome = traverse_circular(adjacency, start, visited)
          chromosomes << ChromosomeSpec.new(blocks: blocks_in_chromosome,
                                            circular: true)
        end

        new(chromosomes,
            blocks: blocks,
            adjacencies: edges.sort_by(&:inspect))
      end

      def initialize(chromosomes, blocks: nil, adjacencies: nil)
        unless chromosomes.is_a?(Array)
          raise ArgumentError, 'a genome must be an array of chromosomes'
        end

        @chromosomes = chromosomes.map { |chromosome| DCJ.coerce_chromosome(chromosome) }
        @blocks = blocks || collect_blocks
        validate!
        @adjacencies = adjacencies || build_adjacencies
      end

      def each(&block)
        @chromosomes.each(&block)
      end

      def size
        @chromosomes.size
      end

      def linear
        @chromosomes.reject(&:circular)
      end

      def circular
        @chromosomes.select(&:circular)
      end

      def adjacencies
        @adjacencies.map(&:dup)
      end

      def block_ids
        @blocks
      end

      def to_a
        @chromosomes.map do |chromosome|
          blocks = chromosome.blocks.map(&DCJ.method(:public_block))
          chromosome.circular ? { blocks: blocks, circular: true } : blocks
        end
      end

      def apply(operation)
        DCJ.apply(self, operation)
      end

      def distance_to(target)
        DCJ.distance(self, target)
      end

      def shortest_sorting_operations(target)
        DCJ.shortest_sorting_operations(self, target)
      end
      alias shortest_path shortest_sorting_operations

      def ==(other)
        other = DCJ.coerce_genome(other)
        normalized_adjacencies == other.send(:normalized_adjacencies)
      rescue ArgumentError
        false
      end
      alias eql? ==

      def hash
        normalized_adjacencies.hash
      end

      private

      def collect_blocks
        @chromosomes.flat_map { |chromosome| chromosome.blocks.map(&:id) }.uniq
      end

      def validate!
        if @chromosomes.empty?
          raise ArgumentError, 'a genome must contain at least one chromosome'
        end

        seen = Set.new
        @chromosomes.each do |chromosome|
          if chromosome.blocks.empty?
            raise ArgumentError, 'chromosomes cannot be empty'
          end

          chromosome.blocks.each do |block|
            unless seen.add?(block.id)
              raise ArgumentError, "duplicate block #{block.id.inspect}"
            end
          end
        end
      end

      def build_adjacencies
        edges = []
        @chromosomes.each do |chromosome|
          blocks = chromosome.blocks
          edges.concat(internal_edges(blocks))
          if chromosome.circular
            edges << edge(end_of(blocks.last, 1), end_of(blocks.first, 0))
          else
            edges << [end_of(blocks.first, 0)]
            edges << [end_of(blocks.last, 1)]
          end
        end
        edges.map do |current_edge|
          self.class.send(:normalize_edge, current_edge)
        end.sort_by(&:inspect)
      end

      def internal_edges(blocks)
        blocks.each_cons(2).map do |left, right|
          edge(end_of(left, 1), end_of(right, 0))
        end
      end

      def edge(left, right)
        [left, right]
      end

      def end_of(block, side)
        if block.orientation == 1
          side.zero? ? [block.id, TAIL] : [block.id, HEAD]
        else
          side.zero? ? [block.id, HEAD] : [block.id, TAIL]
        end
      end

      def normalized_adjacencies
        @adjacencies
      end

      def self.normalize_edge(edge)
        normalized = Array(edge)
        unless normalized.size.between?(1, 2) &&
               normalized.all? { |endpoint| valid_endpoint?(endpoint) }
          raise ArgumentError, "invalid adjacency #{edge.inspect}"
        end
        normalized.sort_by(&:inspect)
      end

      def self.valid_endpoint?(endpoint)
        endpoint.is_a?(Array) && endpoint.size == 2 &&
          [HEAD, TAIL].include?(endpoint.last)
      end

      def self.traverse_linear(adjacency, start, visited)
        blocks = []
        endpoint = start

        loop do
          visited << endpoint
          block_id, side = endpoint
          blocks << SignedBlock.new(id: block_id,
                                   orientation: side == TAIL ? 1 : -1)
          paired_endpoint = [block_id, side == TAIL ? HEAD : TAIL]
          visited << paired_endpoint
          next_edge = adjacency[paired_endpoint]

          break if next_edge.size == 1

          endpoint = other_endpoint(next_edge, paired_endpoint)
        end

        blocks
      end

      def self.traverse_circular(adjacency, start, visited)
        blocks = []
        endpoint = start

        loop do
          visited << endpoint
          block_id, side = endpoint
          blocks << SignedBlock.new(id: block_id,
                                   orientation: side == TAIL ? 1 : -1)
          paired_endpoint = [block_id, side == TAIL ? HEAD : TAIL]
          visited << paired_endpoint
          current_edge = adjacency[paired_endpoint]
          endpoint = other_endpoint(current_edge, paired_endpoint)
          break if endpoint == start
        end

        blocks
      end

      def self.other_endpoint(edge, endpoint)
        raise ArgumentError, 'invalid adjacency set' unless edge

        edge.reject { |candidate| candidate == endpoint }.first
      end

      private_class_method :normalize_edge, :traverse_linear, :traverse_circular,
                           :other_endpoint
    end

    # One valid double-cut-and-join operation.
    class Operation
      attr_reader :cuts, :joins

      def initialize(cuts, joins)
        @cuts = cuts.map { |edge| DCJ::Genome.send(:normalize_edge, edge) }.uniq
        @joins = joins.map { |edge| DCJ::Genome.send(:normalize_edge, edge) }.uniq
        validate!
      end

      def cut
        @cuts
      end

      def join
        @joins
      end

      def to_h
        { cuts: @cuts.map(&:dup), joins: @joins.map(&:dup) }
      end

      def ==(other)
        return false unless other.is_a?(Operation)

        cuts == other.cuts && joins == other.joins
      end
      alias eql? ==

      def hash
        [cuts, joins].hash
      end

      private

      def validate!
        if cuts.empty? || cuts.size > 2 || joins.empty? || joins.size > 2
          raise ArgumentError, 'a DCJ operation must contain one or two cuts and joins'
        end
        if cuts == joins
          raise ArgumentError, 'a DCJ operation must change at least one adjacency'
        end
        unless endpoints(cuts).sort_by(&:inspect) ==
               endpoints(joins).sort_by(&:inspect)
          raise ArgumentError, 'cuts and joins must reconnect the same endpoints'
        end
        unless endpoints(cuts).size == endpoints(cuts).uniq.size
          raise ArgumentError, 'an endpoint cannot be used more than once in an operation'
        end
        unless valid_dcj_shape?
          raise ArgumentError, 'invalid double-cut-and-join operation shape'
        end
      end

      def endpoints(edges)
        edges.flatten(1)
      end

      def valid_dcj_shape?
        cut_internal = cuts.count { |edge| edge.size == 2 }
        cut_telomeres = cuts.count(&:one?)
        join_internal = joins.count { |edge| edge.size == 2 }
        join_telomeres = joins.count(&:one?)
        patterns = [
          [[1, 0], [0, 2]],
          [[0, 2], [1, 0]],
          [[1, 1], [1, 1]],
          [[2, 0], [2, 0]]
        ]
        patterns.include?([[cut_internal, cut_telomeres],
                           [join_internal, join_telomeres]])
      end
    end

    class << self
      def [](*chromosomes)
        Genome.new(chromosomes)
      end

      def genome(chromosomes)
        coerce_genome(chromosomes)
      end

      def operation(cuts, joins)
        Operation.new(cuts, joins)
      end

      def distance(source, target)
        source = coerce_genome(source)
        target = coerce_genome(target)
        compare_blocks!(source, target)

        adjacency_graph = build_adjacency_graph(source, target)
        cycles = 0
        opposite_telomere_paths = 0
        visited = Set.new

        source.blocks.sort_by(&:inspect).each do |block|
          [TAIL, HEAD].each do |side|
            start = [block, side]
            next if visited.include?(start)

            component_endpoints = traverse_component(adjacency_graph, start, visited)
            adjacency_nodes = component_endpoints.select do |node|
              node.is_a?(Array) && [:source, :target].include?(node.first)
            end
            telomere_nodes = adjacency_nodes.select { |node| node[1].size == 1 }
            if telomere_nodes.empty?
              cycles += 1
            elsif telomere_nodes.map(&:first).uniq.size == 2
              opposite_telomere_paths += 1
            end
          end
        end

        source.blocks.size - cycles - opposite_telomere_paths / 2
      end

      def apply(genome, operation)
        genome = coerce_genome(genome)
        unless operation.is_a?(Operation)
          operation = Operation.new(operation[:cuts], operation[:joins])
        end

        current = genome.adjacencies
        operation.cuts.each do |cut|
          unless current.include?(cut)
            raise ArgumentError, "cut #{cut.inspect} is not present in the current genome"
          end
        end

        remaining = current - operation.cuts
        operation.joins.each do |join|
          already_adjacent = remaining.any? { |edge| (edge & join).any? }
          if already_adjacent
            raise ArgumentError, "joined endpoint in #{join.inspect} is already adjacent"
          end
        end

        Genome.from_adjacencies(remaining + operation.joins)
      end

      def shortest_sorting_operations(source, target)
        source = coerce_genome(source)
        target = coerce_genome(target)
        compare_blocks!(source, target)

        operations = []
        current = source

        target.adjacencies.each do |target_edge|
          next if target_edge.size == 1
          next if current.adjacencies.include?(target_edge)

          source_edges = target_edge.map do |endpoint|
            current.adjacencies.find { |edge| edge.include?(endpoint) }
          end
          other_endpoints = source_edges.zip(target_edge).flat_map do |edge, endpoint|
            edge.reject { |candidate| candidate == endpoint }
          end
          joins = [target_edge] + other_endpoints.each_slice(2).to_a
          operation = Operation.new(source_edges.uniq, joins)

          operations << operation
          current = apply(current, operation)
        end

        target_telomere_endpoints = target.adjacencies.select(&:one?).flatten(1)

        current.adjacencies.each do |edge|
          next unless edge.size == 2
          next unless edge.all? { |endpoint| target_telomere_endpoints.include?(endpoint) }

          operation = Operation.new([edge], [[edge.first], [edge.last]])
          operations << operation
          current = apply(current, operation)
        end

        unless current == target
          raise RuntimeError, 'failed to produce a shortest DCJ sorting sequence'
        end
        unless operations.size == distance(source, target)
          raise RuntimeError, 'non-optimal DCJ sorting sequence'
        end

        operations
      end

      def coerce_genome(value)
        return value if value.is_a?(Genome)

        Genome.new(value)
      end

      def coerce_chromosome(value)
        case value
        when ChromosomeSpec
          value
        when Hash
          blocks = value.fetch(:blocks) { value.fetch('blocks') }
          circular = value.fetch(:circular, value.fetch('circular', false))
          ChromosomeSpec.new(blocks: parse_blocks(blocks),
                             circular: circular == true)
        when Array
          ChromosomeSpec.new(blocks: parse_blocks(value), circular: false)
        else
          raise ArgumentError, "invalid chromosome #{value.inspect}"
        end
      end

      def public_block(block)
        return block.id if block.orientation.positive? && block.id.is_a?(Integer)
        return -block.id if block.id.is_a?(Integer)
        "#{block.orientation.positive? ? '+' : '-'}#{block.id}"
      end

      private

      def parse_blocks(value)
        raise ArgumentError, 'blocks must be an array' unless value.is_a?(Array)

        value.map do |block|
          case block
          when SignedBlock
            block
          when Integer
            if block.zero?
              raise ArgumentError, 'block orientation cannot be zero'
            end

            SignedBlock.new(id: block.abs, orientation: block.positive? ? 1 : -1)
          when String
            match = /\A([+-])(.+)\z/.match(block)
            unless match
              raise ArgumentError, "invalid block orientation #{block.inspect}"
            end

            label = /\A\d+\z/.match?(match[2]) ? match[2].to_i : match[2]
            SignedBlock.new(id: label, orientation: match[1] == '+' ? 1 : -1)
          else
            raise ArgumentError, "invalid block #{block.inspect}"
          end
        end
      end

      def compare_blocks!(source, target)
        source_blocks = source.blocks.sort_by(&:inspect)
        target_blocks = target.blocks.sort_by(&:inspect)
        return if source_blocks == target_blocks

        raise ArgumentError, 'genomes must contain the same set of non-repeated blocks'
      end

      def build_adjacency_graph(source, target)
        graph = Hash.new { |hash, key| hash[key] = [] }
        add_adjacency_edges(graph, source.adjacencies, :source)
        add_adjacency_edges(graph, target.adjacencies, :target)
        graph
      end

      def add_adjacency_edges(graph, edges, origin)
        edges.each do |edge|
          adjacency_node = [origin, edge]
          edge.each do |endpoint|
            graph[adjacency_node] << endpoint
            graph[endpoint] << adjacency_node
          end
        end
      end

      def traverse_component(graph, start, visited)
        endpoints = []
        stack = [start]

        until stack.empty?
          node = stack.pop
          next if visited.include?(node)

          visited << node
          endpoints << node
          stack.concat(graph[node])
        end

        endpoints
      end
    end
  end

  # Preferred long name for {DCJ}.
  DoubleCutAndJoin = DCJ
end
