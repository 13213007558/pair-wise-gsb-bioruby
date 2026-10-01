# frozen_string_literal: true

#
# = bio/alignment/graph.rb - alignment of a sequence to a DAG of characters
#
# Copyright:: Copyright (C) 2026 BioRuby project
# License::   The Ruby License
#
# = About Bio::Alignment::Graph and Bio::Alignment::GraphAligner
#
# A set of known variants of a region can be represented as a directed
# acyclic graph (DAG) whose nodes are labeled with single characters,
# instead of collapsing the variants into a single reference sequence.
# Bio::Alignment::Graph holds such a graph (nodes may be shared between
# branches), and Bio::Alignment::GraphAligner aligns a query sequence to
# the best source-to-sink path of the graph with affine gap costs, in
# either global or local mode.
#
# = Usage
#
#   graph = Bio::Alignment::Graph.new
#   a = graph.add_node('A')
#   b = graph.add_node('C')
#   c = graph.add_node('G')
#   d = graph.add_node('T')
#   graph.add_edge(a, b)
#   graph.add_edge(a, c)
#   graph.add_edge(b, d)
#   graph.add_edge(c, d)
#
#   aligner = Bio::Alignment::GraphAligner.new(match: 2, mismatch: -1,
#                                              gap_open: -3, gap_extend: -1,
#                                              mode: :global)
#   result = aligner.align('AGT', graph)
#   result.score            #=> best score over the paths "ACT" and "AGT"
#   result.nodes            #=> graph nodes visited by the traceback
#   result.pairs            #=> [[node_or_nil, sequence_index_or_nil], ...]
#

module Bio
  module Alignment
    # Bio::Alignment::Graph is a directed acyclic graph whose nodes are
    # labeled with single characters. Nodes may be shared between
    # branches (a node may have multiple predecessors and multiple
    # successors). Cycles are rejected: every operation that traverses
    # the graph raises Bio::Alignment::Graph::CyclicError for cyclic
    # input.
    class Graph
      # Raised when the graph contains a cycle.
      class CyclicError < ArgumentError; end

      # Bio::Alignment::Graph::Node is a vertex of Bio::Alignment::Graph
      # labeled with a single character.
      class Node
        # the single-character label of the node
        attr_reader :char
        # position of the node in the insertion order of the graph
        attr_reader :index
        # nodes with an edge pointing to this node
        attr_reader :predecessors
        # nodes this node points to
        attr_reader :successors

        def initialize(char, index) # :nodoc:
          @char = char
          @index = index
          @predecessors = []
          @successors = []
        end

        def to_s
          "#{char}(#{index})"
        end

        def inspect
          "#<#{self.class.name} #{char} index=#{index}>"
        end
      end

      # all nodes of the graph, in insertion order
      attr_reader :nodes

      def initialize
        @nodes = []
      end

      # Adds a new node labeled with the single character +char+
      # and returns it. Raises ArgumentError unless +char+ is exactly
      # one character.
      def add_node(char)
        char = char.to_s
        unless char.length == 1
          raise ArgumentError, "node label must be a single character: #{char.inspect}"
        end

        node = Node.new(char, @nodes.length)
        @nodes << node
        node
      end

      # Adds a directed edge from node +from+ to node +to+.
      # Both nodes must belong to this graph. Returns self.
      def add_edge(from, to)
        check_node(from)
        check_node(to)
        from.successors << to
        to.predecessors << from
        self
      end

      # Creates a chain of new nodes, one per character of +string+,
      # connected by edges in order. Returns the created nodes.
      def add_path(string)
        nodes = string.to_s.each_char.map { |char| add_node(char) }
        nodes.each_cons(2) { |pred, succ| add_edge(pred, succ) }
        nodes
      end

      # nodes without predecessors
      def sources
        @nodes.select { |node| node.predecessors.empty? }
      end

      # nodes without successors
      def sinks
        @nodes.select { |node| node.successors.empty? }
      end

      # Returns the nodes in a topological order.
      # Raises Bio::Alignment::Graph::CyclicError if the graph
      # contains a cycle.
      def topological_order
        order = kahn_sort
        raise CyclicError, 'the graph contains a cycle' if order.size < @nodes.size

        order
      end

      # Yields every source-to-sink path as an array of nodes.
      # Returns an Enumerator when no block is given.
      # Raises Bio::Alignment::Graph::CyclicError for cyclic graphs.
      def each_path(&block)
        return enum_for(:each_path) unless block

        topological_order
        sources.each { |source| walk_path(source, [], &block) }
        self
      end

      # Returns every source-to-sink path as a string of node labels.
      def path_strings
        each_path.map { |path| path.map(&:char).join }
      end

      private

      def check_node(node)
        return if node.is_a?(Node) && @nodes.include?(node)

        raise ArgumentError, 'node does not belong to this graph'
      end

      def kahn_sort
        indegree = @nodes.to_h { |node| [node, node.predecessors.size] }
        queue = @nodes.select { |node| indegree[node].zero? }
        order = []
        until queue.empty?
          node = queue.shift
          order << node
          relax_edges(node, indegree, queue)
        end
        order
      end

      def relax_edges(node, indegree, queue)
        node.successors.each do |succ|
          indegree[succ] -= 1
          queue << succ if indegree[succ].zero?
        end
      end

      def walk_path(node, trail, &block)
        trail = trail + [node]
        return yield(trail) if node.successors.empty?

        node.successors.each { |succ| walk_path(succ, trail, &block) }
      end
    end

    # Bio::Alignment::GraphAlignment is the result of aligning a
    # sequence to a Bio::Alignment::Graph. It reports the score, the
    # graph nodes visited by the traceback, and the corresponding
    # sequence coordinates.
    class GraphAlignment
      # the graph the sequence was aligned to
      attr_reader :graph
      # the aligned query sequence (String)
      attr_reader :sequence
      # alignment mode (:global or :local)
      attr_reader :mode
      # alignment score
      attr_reader :score
      # Alignment columns as [node, sequence_index] pairs, in order.
      # A nil node means a sequence character inserted against a gap in
      # the graph; a nil sequence index means a graph node deleted
      # against a gap in the sequence. Sequence indices are 0-based.
      attr_reader :pairs

      def initialize(graph, sequence, mode, score, pairs) # :nodoc:
        @graph = graph
        @sequence = sequence
        @mode = mode
        @score = score
        @pairs = pairs
      end

      # graph nodes visited by the traceback, in order
      def nodes
        @pairs.map(&:first).compact
      end

      # 0-based index of the first aligned sequence character,
      # or nil when no sequence character is aligned
      def sequence_begin
        sequence_indices.first
      end

      # 0-based end-exclusive index just after the last aligned
      # sequence character, or nil when no sequence character is aligned
      def sequence_end
        index = sequence_indices.last
        index ? index + 1 : nil
      end

      # aligned graph path as a string ('-' for inserted columns)
      def aligned_graph
        @pairs.map { |node, _| node ? node.char : '-' }.join
      end

      # aligned query sequence as a string ('-' for deleted columns)
      def aligned_sequence
        @pairs.map { |_, index| index ? @sequence[index] : '-' }.join
      end

      def to_s
        "#{aligned_graph}\n#{aligned_sequence}"
      end

      private

      def sequence_indices
        @pairs.map(&:last).compact
      end
    end

    # Bio::Alignment::GraphAligner aligns a query sequence to a
    # Bio::Alignment::Graph (a DAG of character-labeled nodes) by
    # dynamic programming over the topological order of the graph.
    #
    # Both global (Needleman-Wunsch-like: the whole sequence against a
    # complete source-to-sink path) and local (Smith-Waterman-like)
    # modes are supported, with a substitution score and affine gap
    # costs. A gap of length L costs gap_open + (L - 1) * gap_extend,
    # so a continuous gap spanning several graph nodes is opened only
    # once.
    #
    # Tie-breaking is deterministic: when several predecessors or
    # states give the same score, the predecessor node inserted into
    # the graph first wins, match/mismatch states are preferred over
    # gap states, and extending a gap is preferred over opening a new
    # one. Parallel paths with equal scores therefore always yield the
    # same traceback.
    class GraphAligner
      # score for a pair of equal characters
      attr_reader :match
      # score for a pair of different characters
      attr_reader :mismatch
      # optional substitution table (a Hash keyed by two-character
      # arrays, or a callable taking two characters)
      attr_reader :substitution
      # cost of opening a gap
      attr_reader :gap_open
      # cost of extending a gap by one character
      attr_reader :gap_extend
      # alignment mode (:global or :local)
      attr_reader :mode

      def initialize(match: 2, mismatch: -1, substitution: nil,
                     gap_open: -3, gap_extend: -1, mode: :global)
        @match = match
        @mismatch = mismatch
        @substitution = substitution
        @gap_open = gap_open
        @gap_extend = gap_extend
        @mode = mode.to_sym
        return if %i[global local].include?(@mode)

        raise ArgumentError, "invalid mode: #{mode.inspect} (must be :global or :local)"
      end

      # Aligns +sequence+ to +graph+ and returns a
      # Bio::Alignment::GraphAlignment. Raises
      # Bio::Alignment::Graph::CyclicError when the graph has a cycle,
      # and ArgumentError when the graph has no nodes.
      def align(sequence, graph)
        seq = sequence.to_s
        raise ArgumentError, 'the graph has no nodes' if graph.nodes.empty?

        GraphAlignerEngine.new(self, graph, seq).run
      end

      # substitution score between two characters
      def substitution_score(char_a, char_b)
        return char_a == char_b ? @match : @mismatch unless @substitution

        lookup_substitution(char_a, char_b)
      end

      private

      def lookup_substitution(char_a, char_b)
        return @substitution.call(char_a, char_b) if @substitution.respond_to?(:call)

        @substitution.fetch([char_a, char_b]) do
          @substitution.fetch([char_b, char_a]) do
            char_a == char_b ? @match : @mismatch
          end
        end
      end
    end

    # Bio::Alignment::GraphAlignerEngine is the dynamic-programming
    # engine used internally by Bio::Alignment::GraphAligner.
    # It is not part of the public API.
    class GraphAlignerEngine # :nodoc:
      STATES = %i[M Ix Iy].freeze
      NEGATIVE_INFINITY = -Float::INFINITY
      # pseudo-node preceding all source nodes
      START = Object.new

      def initialize(aligner, graph, sequence)
        @aligner = aligner
        @graph = graph
        @seq = sequence
        @len = sequence.length
        @scores = {}
        @pointers = {}
      end

      def run
        init_start_rows
        @graph.topological_order.each { |node| compute_node(node) }
        node, state, pos = find_end
        GraphAlignment.new(@graph, @seq, @aligner.mode, @end_score,
                           traceback(node, state, pos))
      end

      private

      def init_start_rows
        if @aligner.mode == :global
          set_row(START, :M, [0] + [NEGATIVE_INFINITY] * @len)
          set_row(START, :Ix, insertion_row)
          set_row(START, :Iy, [NEGATIVE_INFINITY] * (@len + 1))
          (1..@len).each { |j| @pointers[[START, :Ix]][j] = [START, j == 1 ? :M : :Ix] }
        else
          STATES.each { |state| set_row(START, state, [0] * (@len + 1)) }
        end
      end

      def insertion_row
        (0..@len).map do |j|
          j.zero? ? NEGATIVE_INFINITY : @aligner.gap_open + (j - 1) * @aligner.gap_extend
        end
      end

      def set_row(node, state, row)
        @scores[[node, state]] = row
        @pointers[[node, state]] = Array.new(@len + 1)
      end

      def compute_node(node)
        preds = predecessors_of(node)
        STATES.each { |state| set_row(node, state, Array.new(@len + 1)) }
        (0..@len).each { |j| compute_cell(node, preds, j) }
      end

      def compute_cell(node, preds, j)
        store(node, :Ix, j, choose(ix_candidates(node, j))) if j.positive?
        store(node, :Iy, j, choose(iy_candidates(preds, j)))
        store(node, :M, j, choose(m_candidates(node, preds, j))) if j.positive?
      end

      def m_candidates(node, preds, j)
        preds.flat_map do |pred|
          STATES.map do |state|
            [score_at(pred, state, j - 1) + substitution(node, j), [pred, state]]
          end
        end
      end

      def ix_candidates(node, j)
        [[score_at(node, :Ix, j - 1) + @aligner.gap_extend, [node, :Ix]],
         [score_at(node, :M, j - 1) + @aligner.gap_open, [node, :M]]]
      end

      def iy_candidates(preds, j)
        preds.map { |p| [score_at(p, :Iy, j) + @aligner.gap_extend, [p, :Iy]] } +
          preds.map { |p| [score_at(p, :M, j) + @aligner.gap_open, [p, :M]] }
      end

      def substitution(node, j)
        @aligner.substitution_score(node.char, @seq[j - 1])
      end

      def predecessors_of(node)
        preds = node.predecessors.sort_by(&:index)
        preds.empty? ? [START] : preds
      end

      def choose(candidates)
        best = candidates.first
        candidates.each { |cand| best = cand if cand[0] > best[0] }
        return best unless @aligner.mode == :local

        best[0].positive? ? best : [0, nil]
      end

      def store(node, state, j, best)
        @scores[[node, state]][j] = best[0]
        @pointers[[node, state]][j] = best[1]
      end

      def score_at(node, state, j)
        row = @scores[[node, state]]
        value = row && row[j]
        value || NEGATIVE_INFINITY
      end

      def find_end
        @aligner.mode == :global ? find_end_global : find_end_local
      end

      def find_end_global
        candidates = @graph.sinks.sort_by(&:index).flat_map do |node|
          STATES.map { |state| [score_at(node, state, @len), [node, state]] }
        end
        best = candidates.first
        candidates.each { |cand| best = cand if cand[0] > best[0] }
        @end_score = best[0]
        [best[1][0], best[1][1], @len]
      end

      def find_end_local
        best = [0, nil, :M, 0]
        @graph.nodes.sort_by(&:index).each do |node|
          best = local_best_at(node, best)
        end
        @end_score = best[0]
        best[1..]
      end

      def local_best_at(node, best)
        (0..@len).each do |j|
          STATES.each do |state|
            score = score_at(node, state, j)
            best = [score, node, state, j] if score > best[0]
          end
        end
        best
      end

      def traceback(node, state, j)
        pairs = []
        while node
          break if @aligner.mode == :local && score_at(node, state, j) <= 0

          pair = emit_pair(node, state, j)
          pairs << pair if pair
          node, state, j = step_back(node, state, j)
        end
        pairs.reverse
      end

      def emit_pair(node, state, j)
        case state
        when :Ix then [nil, j - 1]
        when :M then node.equal?(START) ? nil : [node, j - 1]
        else [node, nil]
        end
      end

      def step_back(node, state, j)
        prev = @pointers[[node, state]][j]
        return [nil, nil, nil] unless prev

        [prev[0], prev[1], state == :Iy ? j : j - 1]
      end
    end
  end
end
