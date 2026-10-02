# frozen_string_literal: true

module Bio
  # Nucleic-acid assembly utilities implemented without external programs.
  module Assembly
    # Canonical bidirectional de Bruijn graph for short nucleotide reads.
    class DeBruijnGraph
      attr_reader :k, :kmer_counts, :edge_counts

      def self.build(reads, k)
        graph = new(k)
        reads.each { |read| graph.add_read(read) }
        graph
      end

      def self.from_reads(reads, k)
        build(reads, k)
      end

      def initialize(k, reads = nil)
        raise ArgumentError, 'k must be a positive integer' unless
          k.is_a?(Integer) && k.positive?

        @k = k
        @kmer_counts = Hash.new(0)
        @edge_counts = Hash.new(0)
        @outgoing = Hash.new { |hash, side| hash[side] = {} }
        @incoming = Hash.new { |hash, side| hash[side] = {} }
        @sides = {}

        reads.to_a.each { |read| add_read(read) } if reads
      end

      def add_read(read)
        sequence = read.to_s.downcase.tr('u', 't').gsub(/\s+/, '')
        sequence.split(/[^acgt]/).each { |window| add_window(window) }
        self
      end

      def unitigs
        paths, processed = linear_unitigs
        paths.concat(remaining_unitigs(processed))
        paths.uniq(&:unitig_key).sort_by do |unitig|
          [
            unitig.sequence,
            unitig.start_node,
            unitig.start_orientation,
            unitig.end_node,
            unitig.end_orientation
          ]
        end
      end

      def outgoing_states(side)
        @outgoing[side].keys.sort
      end

      def incoming_states(side)
        @incoming[side].keys.sort
      end

      def canonical_kmer(kmer)
        complement = reverse_complement(kmer)
        [kmer, complement].min
      end

      private

      def add_window(window)
        return if window.length < @k

        0.upto(window.length - @k) do |position|
          kmer = window[position, @k]
          node = canonical_kmer(kmer)
          @kmer_counts[node] += 1
          register_node(node)

          next if position + @k >= window.length

          add_edge(window[position, @k + 1])
        end
      end

      def register_node(node)
        @sides[[node, '+']] = true
        @sides[[node, '-']] = true
      end

      def add_edge(edge)
        edge_id = canonical_kmer(edge)

        prefix = edge[0, @k]
        suffix = edge[1, @k]
        from = side_for_kmer(prefix)
        to = side_for_kmer(suffix)
        reverse_from = opposite_side(to)
        reverse_to = opposite_side(from)
        add_arc(from, to, edge_id)
        @edge_counts[edge_id] += 1

        add_arc(reverse_from, reverse_to, edge_id)
      end

      def add_arc(from, to, edge)
        @outgoing[from][to] = edge
        @incoming[to][from] = edge
      end

      def side_for_kmer(kmer)
        node = canonical_kmer(kmer)
        [node, kmer == node ? '+' : '-']
      end

      def opposite_side(side)
        node, orientation = side
        [node, orientation == '+' ? '-' : '+']
      end

      def side_kmer(side)
        node, orientation = side
        orientation == '+' ? node : reverse_complement(node)
      end

      def sorted_sides
        @sides.keys.sort
      end

      def outgoing_edges(side)
        @outgoing[side].values.uniq
      end

      def incoming_edges(side)
        @incoming[side].values.uniq
      end

      def middle_side?(side)
        incoming_edges(side).size == 1 && outgoing_edges(side).size == 1
      end

      def linear_unitigs
        result = []
        processed = {}

        sorted_sides.each do |side|
          next if middle_side?(side)

          outgoing_states(side).each do |target|
            path = walk(side, target)
            path[:states].each { |item| processed[item] = true }
            result << build_unitig(
              path[:states], path[:edges], false
            ) if path[:edges].any?
          end
        end

        singleton_unitigs(processed, result)
        [result, processed]
      end

      def walk(start, first_target)
        states = [start]
        edges = []
        used_edges = {}
        visited_sides = { start => true }
        current = start
        target = first_target

        loop do
          edge = @outgoing[current][target]
          break if edge.nil? || used_edges.key?(edge)

          used_edges[edge] = true
          edges << edge
          states << target
          current = target
          break if current == start
          break unless middle_side?(current)
          break if visited_sides.key?(current)

          visited_sides[current] = true
          target = outgoing_states(current).first
        end

        { states: states, edges: edges }
      end

      def singleton_unitigs(processed, result)
        emitted = {}

        sorted_sides.each do |node, orientation|
          side = [node, orientation]
          next unless @outgoing[side].empty? && @incoming[side].empty?
          next if emitted.key?(node) || processed.key?(side)

          other_side = [node, orientation == '+' ? '-' : '+']
          next unless @outgoing[other_side].empty? &&
                      @incoming[other_side].empty?

          emitted[node] = true
          result << build_unitig([[node, '+']], [], false)
        end
      end

      def remaining_unitigs(processed)
        result = []

        sorted_sides.each do |side|
          next unless middle_side?(side)
          next if processed.key?(side)

          path = remaining_walk(side)
          path[:states].each { |item| processed[item] = true }
          path[:states].each { |item| processed[opposite_side(item)] = true } if
            path[:circular]
          next if path[:edges].empty?

          result << build_unitig(
            path[:states], path[:edges], path[:circular]
          )
        end

        result
      end

      def remaining_walk(start)
        states = [start]
        edges = []
        used_edges = {}
        positions = { start => 0 }
        current = start

        while middle_side?(current)
          target = outgoing_states(current).first
          edge = @outgoing[current][target]
          return { states: states, edges: edges, circular: false } if
            used_edges.key?(edge)
          if target == start
            edges << edge
            states << target
            return { states: states, edges: edges, circular: true }
          end
          return { states: states, edges: edges, circular: false } if
            positions.key?(target)

          used_edges[edge] = true
          edges << edge
          states << target
          positions[target] = states.size - 1 unless target == start
          current = target
        end

        { states: states, edges: edges, circular: false }
      end

      def build_unitig(states, edges, circular)
        forward_sequence = spell(states, false)
        reverse_states = states.reverse.map { |side| opposite_side(side) }
        reverse_sequence = spell(reverse_states, false)
        forward_key = states.map do |node, orientation|
          orientation == '-' ? reverse_complement(node) : node
        end
        reverse_key = reverse_states.map do |node, orientation|
          orientation == '-' ? reverse_complement(node) : node
        end

        if circular
          states, edges = canonical_circular_rotation(states, edges)
          sequence = spell(states, true)
        else
          orientation = [reverse_sequence, reverse_key] <=>
                        [forward_sequence, forward_key]
          if orientation == -1
            states = reverse_states
            edges = edges.reverse
            sequence = reverse_sequence
          else
            sequence = forward_sequence
          end
        end

        Unitig.new(sequence, states, edges.map { |edge| @edge_counts[edge] },
                   @kmer_counts, @k, circular)
      end

      def spell(states, circular)
        sequence = side_kmer(states.first).dup
        trailing_states = circular ? states[1...-1] : states[1..]
        trailing_states.each do |side|
          sequence << side_kmer(side)[-1]
        end
        sequence
      end

      def canonical_circular_rotation(states, edges)
        cycle_states = states[0...-1]
        candidates = circular_candidates(cycle_states, edges, 0)
        reversed_states = cycle_states.reverse.map do |side|
          opposite_side(side)
        end
        candidates.concat(
          circular_candidates(reversed_states, edges.reverse, 1)
        )
        candidate = candidates.min_by { |item| item[0, 3] }
        [candidate[3] + [candidate[3].first], candidate[4]]
      end

      def circular_candidates(states, edges, direction)
        states.size.times.map do |rotation|
          rotated = states.rotate(rotation)
          closure = rotated + [rotated.first]
          [spell(closure, true), direction, rotation, rotated,
           edges.rotate(rotation)]
        end
      end

      def reverse_complement(sequence)
        sequence.tr('acgt', 'tgca').reverse
      end

      # A maximal non-branching path in a bidirectional de Bruijn graph.
      class Unitig
        attr_reader :sequence, :kmer_counts, :edge_counts, :k,
                    :start_kmer, :end_kmer, :start_node, :end_node,
                    :start_orientation, :end_orientation

        def initialize(sequence, states, edge_counts, kmer_counts, k,
                       circular = false)
          @sequence = sequence.dup
          @states = states.map(&:dup)
          @edge_counts = edge_counts.dup
          @k = k
          @circular = circular
          @start_node, @start_orientation = @states.first
          @end_node, @end_orientation = @states.last
          @start_kmer = state_kmer(@states.first)
          @end_kmer = state_kmer(@states.last)
          coverage_states = circular ? @states[0...-1] : @states
          @kmer_counts = coverage_states.map do |node, _orientation|
            kmer_counts[node]
          end
        end

        def circular?
          @circular
        end

        def mean_kmer_coverage
          mean(@kmer_counts)
        end

        def mean_edge_coverage
          mean(@edge_counts)
        end

        def coverage
          counts = @edge_counts.empty? ? @kmer_counts : @edge_counts
          mean(counts).to_f
        end

        def unitig_key
          end_node_string = canonical_kmer(end_kmer)
          end_sign = end_kmer == end_node_string ? '+' : '-'
          start_rc = reverse_complement(start_kmer)
          start_node_string = canonical_kmer(start_rc)
          end_half_edge = [end_node_string, end_sign]
          start_sign = start_rc == start_node_string ? '+' : '-'
          start_half_edge = [start_node_string, start_sign]
          canonical_endpoint = [end_half_edge, start_half_edge].sort
          canonical_sequence = @sequence <= reverse_complement(@sequence) ?
                               @sequence : reverse_complement(@sequence)
          [
            canonical_sequence, @circular, @kmer_counts, @edge_counts,
            canonical_endpoint
          ]
        end

        def ==(other)
          other.is_a?(Unitig) && unitig_key == other.unitig_key
        end
        alias eql? ==

        def hash
          unitig_key.hash
        end

        def to_s
          @sequence
        end

        private

        def mean(counts)
          return 0 if counts.empty?

          Rational(counts.inject(0, :+), counts.size)
        end

        def state_kmer(side)
          node, orientation = side
          orientation == '+' ? node : reverse_complement(node)
        end

        def canonical_kmer(kmer)
          complement = reverse_complement(kmer)
          [kmer, complement].min
        end

        def reverse_complement(sequence)
          sequence.tr('acgt', 'tgca').reverse
        end
      end
    end

    def self.unitigs(reads, k)
      DeBruijnGraph.build(reads, k).unitigs
    end

    def self.build_unitigs(reads, k)
      unitigs(reads, k)
    end
  end
end
