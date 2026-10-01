# frozen_string_literal: true

#
# = bio/dag_alignment.rb - Sequence to directed acyclic character graph alignment
#
# Copyright:: Copyright (C) 2026 BioRuby project
# License:: The Ruby License
#

module Bio

# Bio::DAGAlignment aligns a linear sequence against a directed acyclic
# graph (DAG) whose nodes each carry exactly one character.  Graph paths
# represent all known variants of a region, including branches and
# shared (re-converging) nodes, without collapsing them onto one
# reference sequence.
#
# Both global (Needleman-Wunsch style) and local (Smith-Waterman style)
# alignments are supported with a substitution score and affine gap
# costs (a gap of length +k+ costs +gap_open + (k - 1) * gap_extend+,
# where both values are normally negative).  The dynamic program keeps
# three states per (node, sequence position) cell: match/mismatch (M),
# graph character aligned to a gap (Y) and sequence character aligned
# to a gap immediately before a node (X).
#
# A gap that continues across a graph edge (a graph-gap run spanning
# several nodes, or a sequence-gap run emitted while crossing an edge)
# is charged the opening penalty exactly once; only the extension
# penalty is paid for subsequent gap columns.
#
# Edges are the only way to move between nodes during fill and
# traceback: two branches that merely contain identical characters are
# never spliced together.
#
# When several alignments share the optimal score, the result is
# chosen deterministically.  At every tie, candidates are considered in
# the fixed state order M, Y, X, and predecessor nodes in topological
# order (ties resolved by insertion order); the first candidate with
# the best score wins.
module DAGAlignment

  # Default score for two equal characters.
  DEFAULT_MATCH = 2

  # Default score for two different characters.
  DEFAULT_MISMATCH = -1

  # Default cost charged once when a gap opens.
  DEFAULT_GAP_OPEN = -2

  # Default cost charged for every additional column of an open gap.
  DEFAULT_GAP_EXTEND = -1

  # Base error class for Bio::DAGAlignment.
  class Error < StandardError
  end

  # Raised when a graph contains a directed cycle or a self-loop.
  class CyclicGraphError < Error
  end

  # Directed acyclic graph whose nodes each carry exactly one character.
  #
  # Node identities may be any hashable objects; edges are given as
  # [from, to] pairs of node identities.  Nodes may be shared by
  # several branches.  Construction raises CyclicGraphError if the edges
  # contain a directed cycle or a self-loop.
  #
  #   graph = Bio::DAGAlignment::CharacterDAG.new(
  #     { 1 => 'A', 2 => 'C', 3 => 'G', 4 => 'T' },
  #     [[1, 2], [2, 4], [1, 3], [3, 4]]
  #   )
  class CharacterDAG
    # Array of node identities in insertion order.
    attr_reader :nodes
    # Array of single-character labels, parallel to #nodes.
    attr_reader :labels
    # Array of [from, to] node-identity pairs.
    attr_reader :edges

    # Builds a DAG from a node-identity => label map and edge pairs.
    #
    # Every endpoint in +edges+ must be present in +labels+, and every
    # label must be a single-character String.
    def initialize(labels, edges = [])
      unless labels.is_a?(Hash)
        raise ArgumentError, 'labels must be a Hash of node => character'
      end

      @nodes = labels.keys
      @labels = @nodes.map { |node| validate_label(node, labels[node]) }
      @node_index = {}
      @nodes.each_with_index { |node, idx| @node_index[node] = idx }
      @edges = edges.map { |from, to| validate_edge(from, to) }

      @adjacency = Array.new(@nodes.size) { [] }
      @indegree = Array.new(@nodes.size, 0)
      build_adjacency
      @topological_order = topological_sort
    end

    # Number of nodes.
    def size
      @nodes.size
    end

    # Returns the single-character label of node index +index+.
    def label_at(index)
      @labels[index]
    end

    # Returns the node identity stored at internal index +index+.
    def node_at(index)
      @nodes[index]
    end

    # Returns the internal index of node identity +node+, or nil.
    def index_of(node)
      @node_index[node]
    end

    # Returns internal indices in deterministic topological order.
    def topological_order
      @topological_order.dup
    end

    # Returns the successor node indices of +index+.
    def successors(index)
      @adjacency[index].dup
    end

    # Returns the predecessor node indices of +index+.
    def predecessors(index)
      @topological_order.select do |candidate|
        @adjacency[candidate].include?(index)
      end
    end

    # Returns internal indices of nodes without incoming edges.
    def source_indices
      @nodes.each_index.select { |index| @indegree[index].zero? }
    end

    # Returns internal indices of nodes without outgoing edges.
    def sink_indices
      @nodes.each_index.select { |index| @adjacency[index].empty? }
    end

    # Yields every source-to-sink path as an array of node identities.
    # Returns an Enumerator when called without a block.
    def each_path
      return enum_for(:each_path) unless block_given?

      source_indices.each do |source|
        extend_path(source, [source]) do |path|
          yield path.map { |index| @nodes[index] }
        end
      end
      nil
    end

    # Returns the concatenated label string of every source-to-sink path.
    def path_strings
      each_path.map do |path|
        path.map { |node| @labels[@node_index[node]] }.join
      end
    end

    private

    def validate_label(node, value)
      unless value.is_a?(String) && value.length == 1
        raise ArgumentError, "label of node #{node.inspect} must be a single-character String"
      end
      value
    end

    def validate_edge(from, to)
      unless @node_index.key?(from) && @node_index.key?(to)
        raise ArgumentError, "edge #{from.inspect}->#{to.inspect} references an unknown node"
      end
      [from, to]
    end

    def build_adjacency
      @edges.each do |from, to|
        from_idx = @node_index[from]
        to_idx = @node_index[to]
        next if @adjacency[from_idx].include?(to_idx)

        @adjacency[from_idx] << to_idx
        @indegree[to_idx] += 1
      end
    end

    # Kahn topological sort.  Ready nodes are processed in ascending
    # index order so the result is independent of Hash iteration order.
    def topological_sort
      indegree = @indegree.dup
      ready = indegree.each_index.select { |index| indegree[index].zero? }
      order = []
      until ready.empty?
        index = ready.shift
        order << index
        @adjacency[index].sort.each do |successor|
          indegree[successor] -= 1
          ready << successor if indegree[successor].zero?
        end
        ready.sort!
      end
      unless order.size == @nodes.size
        raise CyclicGraphError,
              'graph must be acyclic; a directed cycle (or self-loop) was found'
      end
      order
    end

    def extend_path(index, prefix)
      successors = @adjacency[index]
      if successors.empty?
        yield prefix
      else
        successors.each do |successor|
          extend_path(successor, prefix + [successor]) { |path| yield path }
        end
      end
    end
  end

  # One aligned column produced by traceback.
  #
  # +state+ is +:match+, +:graph_gap+ (graph character versus sequence
  # gap) or +:sequence_gap+ (sequence character versus graph gap).
  # +node+ is the traversed graph node identity, or nil for sequence
  # gaps; +node_index+ and +seq_index+ are zero-based coordinates of
  # the aligned characters, or nil on the gapped side; +score+ is the
  # score contributed by the column.
  class Step
    # Returns the alignment state symbol.
    attr_reader :state
    # Returns the traversed graph node identity, or nil for sequence gaps.
    attr_reader :node
    # Returns the zero-based graph node index, or nil for sequence gaps.
    attr_reader :node_index
    # Returns the zero-based sequence index, or nil for graph gaps.
    attr_reader :seq_index
    # Returns the score contributed by this column.
    attr_reader :score

    # Creates a new traceback step.
    def initialize(state, node, node_index, seq_index, score)
      @state = state
      @node = node
      @node_index = node_index
      @seq_index = seq_index
      @score = score
    end
  end

  # Result of a sequence-to-DAG alignment.
  class Result
    # Total alignment score.
    attr_reader :score
    # Ordered traceback columns (Step objects, forward in alignment order).
    attr_reader :trace
    # Identities of graph nodes visited by the alignment, in order.
    attr_reader :graph_path
    # Path characters with '-' inserted at sequence-gap columns.
    attr_reader :aligned_graph
    # Sequence characters with '-' inserted at graph-gap columns.
    attr_reader :aligned_sequence
    # Zero-based sequence index of the first aligned column (nil if empty).
    attr_reader :start_seq_index
    # Exclusive zero-based sequence index after the last column (nil if empty).
    attr_reader :end_seq_index
    # Identity of the first visited graph node (nil if empty).
    attr_reader :start_node
    # Identity of the last visited graph node (nil if empty).
    attr_reader :end_node

    # Builds a result.  +trace+ must already be in forward order.
    def initialize(score, trace, graph, sequence)
      @score = score
      @trace = trace
      @graph = graph
      @graph_path = trace.filter_map(&:node)
      @aligned_graph = trace.map { |step| graph_column(step) }.join
      @aligned_sequence = trace.map { |step| sequence_column(step, sequence) }.join
      fill_sequence_range
      fill_node_range
    end

    # Returns true when the alignment contains no columns.
    def empty?
      @trace.empty?
    end

    private

    def graph_column(step)
      step.state == :sequence_gap ? '-' : @graph.label_at(step.node_index)
    end

    def sequence_column(step, sequence)
      step.seq_index.nil? ? '-' : sequence[step.seq_index]
    end

    def fill_sequence_range
      matched = @trace.reject { |step| step.seq_index.nil? }
      @start_seq_index = matched.first&.seq_index
      @end_seq_index = matched.empty? ? nil : matched.last.seq_index + 1
    end

    def fill_node_range
      nodes = @trace.reject { |step| step.node_index.nil? }
      @start_node = nodes.first&.node
      @end_node = nodes.last&.node
    end
  end

  # Internal index of the virtual start row stored alongside node rows.
  START_ROW = -1

  # Pointer sentinel for the global alignment start at the virtual row.
  POINTER_START = :start

  # Pointer sentinel for a free local-alignment restart (score 0).
  POINTER_LOCAL = :local

  NEGATIVE_INFINITY = -Float::INFINITY
  private_constant :NEGATIVE_INFINITY

  # Sequence-to-DAG aligner.  Usually invoked through
  # Bio::DAGAlignment.align rather than directly.
  class Aligner
    # Creates an aligner.
    #
    # +mode+ is +:global+ or +:local+.  +match+ and +mismatch+ are the
    # substitution scores for equal/different characters; alternatively
    # a +substitution+ block taking two characters can be given.
    # +gap_open+ is charged for the first column of every gap and
    # +gap_extend+ for each following column; both should be negative.
    def initialize(graph, sequence, mode, match, mismatch, gap_open, gap_extend, &substitution)
      @graph = graph
      @sequence = sequence
      @mode = mode
      @match = match || DEFAULT_MATCH
      @mismatch = mismatch || DEFAULT_MISMATCH
      @gap_open = gap_open
      @gap_extend = gap_extend
      @substitution = substitution
    end

    # Runs the dynamic program and returns a Result.
    def run
      validate
      initialize_tables
      fill
      start_cell = select_start_cell
      return empty_result if start_cell.nil?

      trace = traceback(start_cell)
      Result.new(start_cell[:score], trace, @graph, @sequence)
    end

    private

    def validate
      unless @graph.is_a?(CharacterDAG)
        raise ArgumentError, 'graph must be a Bio::DAGAlignment::CharacterDAG'
      end
      unless @sequence.is_a?(String)
        raise ArgumentError, 'sequence must be a String'
      end
      unless %i[global local].include?(@mode)
        raise ArgumentError, "mode must be :global or :local, not #{@mode.inspect}"
      end
      if @graph.size.zero?
        raise ArgumentError, 'cannot align against an empty graph'
      end
      unless @gap_open <= @gap_extend
        raise ArgumentError, 'gap_open must be <= gap_extend (opening costs at least a continuation)'
      end
    end

    def initialize_tables
      @n = @sequence.length
      @node_count = @graph.size
      rows = @node_count + 1
      @m = Array.new(rows) { Array.new(@n + 1, NEGATIVE_INFINITY) }
      @y = Array.new(rows) { Array.new(@n + 1, NEGATIVE_INFINITY) }
      @x = Array.new(rows) { Array.new(@n + 1, NEGATIVE_INFINITY) }
      @pm = Array.new(rows) { Array.new(@n + 1) }
      @py = Array.new(rows) { Array.new(@n + 1) }
      @px = Array.new(rows) { Array.new(@n + 1) }
      @m[0][0] = 0
      @pm[0][0] = POINTER_START
      @incoming_rows = Array.new(rows) { [] }
      @graph.topological_order.each do |node_index|
        @incoming_rows[row_of(node_index)] = incoming_rows(node_index)
      end
      fill_start_insertions
    end

    # Leading insertions before any graph node, hosted on virtual row 0.
    def fill_start_insertions
      (1..@n).each do |seq_index|
        update(:x, 0, seq_index,
               @x[0][seq_index - 1] + @gap_extend, [:x, 0, seq_index - 1])
        update(:x, 0, seq_index,
               @m[0][seq_index - 1] + @gap_open, [:m, 0, seq_index - 1])
      end
    end

    # Virtual start row is index 0; node v is stored at row v + 1.
    def row_of(node_index)
      node_index + 1
    end

    def fill
      @graph.topological_order.each { |node_index| fill_node(node_index) }
    end

    def fill_node(node_index)
      row = row_of(node_index)
      fill_graph_gap(row, 0)
      (1..@n).each do |seq_index|
        fill_insertions(row, seq_index)
        fill_graph_gap(row, seq_index)
        fill_matches(node_index, row, seq_index)
      end
    end

    # Incoming rows for a node: virtual start row first (sources only),
    # then predecessors in topological order.  This iteration order is
    # what makes equal-score paths resolve deterministically.
    def incoming_rows(node_index)
      rows = []
      rows << 0 if @graph.source_indices.include?(node_index)
      @graph.predecessors(node_index).each { |pred| rows << row_of(pred) }
      rows
    end

    # X: sequence characters aligned to gaps immediately before this
    # node.  The run opens over an incoming edge and continues at this
    # node, so its opening penalty is paid only once for the whole run.
    def fill_insertions(row, seq_index)
      update(:x, row, seq_index,
             @x[row][seq_index - 1] + @gap_extend, [:x, row, seq_index - 1])
      update(:x, row, seq_index,
             @m[row][seq_index - 1] + @gap_open, [:m, row, seq_index - 1])
    end

    # Y: this graph character aligns to a sequence gap.  A deletion run
    # continuing over an edge only pays the extension penalty.
    def fill_graph_gap(row, seq_index)
      @incoming_rows[row].each do |pred_row|
        update(:y, row, seq_index,
               @y[pred_row][seq_index] + @gap_extend, [:y, pred_row, seq_index])
      end
      @incoming_rows[row].each do |pred_row|
        update(:y, row, seq_index,
               @m[pred_row][seq_index] + @gap_open, [:m, pred_row, seq_index])
      end
    end

    # M: this graph character matches (or mismatches) a sequence char.
    def fill_matches(node_index, row, seq_index)
      sub = substitution_score(@graph.label_at(node_index), @sequence[seq_index - 1])
      @incoming_rows[row].each do |pred_row|
        update(:m, row, seq_index,
               @m[pred_row][seq_index - 1] + sub, [:m, pred_row, seq_index - 1])
      end
      update(:m, row, seq_index,
             @x[row][seq_index - 1] + sub, [:x, row, seq_index - 1])
      update(:m, row, seq_index, 0, POINTER_LOCAL) if @mode == :local
    end

    def substitution_score(graph_char, seq_char)
      if @substitution
        @substitution.call(graph_char, seq_char)
      elsif graph_char == seq_char
        @match
      else
        @mismatch
      end
    end

    # Keeps the best (score, pointer) candidate.  Equality keeps the
    # first candidate, which gives deterministic M/Y/X and topological
    # predecessor order tie breaking.
    def update(state, row, seq_index, score, new_pointer)
      return if score == NEGATIVE_INFINITY || score.nan?

      table = state_table(state)
      pointers = pointer_table(state)
      return unless score > table[row][seq_index]

      table[row][seq_index] = score
      pointers[row][seq_index] = new_pointer
    end

    def state_table(state)
      { m: @m, y: @y, x: @x }.fetch(state)
    end

    def pointer_table(state)
      { m: @pm, y: @py, x: @px }.fetch(state)
    end

    def select_start_cell
      best = nil
      candidates.each do |candidate|
        best = candidate if best.nil? || candidate[:score] > best[:score]
      end
      best
    end

    def candidates
      if @mode == :global
        global_candidates
      else
        best_local = best_local_match
        best_local ? [best_local] : []
      end
    end

    def global_candidates
      @graph.sink_indices.flat_map do |node_index|
        row = row_of(node_index)
        [[:m, @m[row][@n]], [:y, @y[row][@n]], [:x, @x[row][@n]]].map do |state, score|
          { state: state, node: node_index, seq: @n, score: score }
        end
      end
    end

    def best_local_match
      best = nil
      @graph.topological_order.each do |node_index|
        row = row_of(node_index)
        (1..@n).each do |seq_index|
          score = @m[row][seq_index]
          next if score.negative?

          candidate = { state: :m, node: node_index, seq: seq_index, score: score }
          best = candidate if best.nil? || score > best[:score]
        end
      end
      best
    end

    def empty_result
      Result.new(0, [], @graph, @sequence)
    end

    def traceback(cell)
      steps = []
      state = cell[:state]
      row = row_of(cell[:node])
      seq_index = cell[:seq]
      loop do
        ptr = pointer_table(state)[row][seq_index]
        break if ptr == POINTER_START || ptr == POINTER_LOCAL || ptr.nil?

        steps << build_step(state, row, seq_index, ptr)
        state, row, seq_index = ptr
      end
      steps.reverse
    end

    def build_step(state, row, seq_index, ptr)
      case state
      when :m
        node_index = row - 1
        score = @m[row][seq_index] - predecessor_score(ptr)
        Step.new(:match, @graph.node_at(node_index), node_index, seq_index - 1, score)
      when :y
        node_index = row - 1
        score = @y[row][seq_index] - predecessor_score(ptr)
        Step.new(:graph_gap, @graph.node_at(node_index), node_index, nil, score)
      else
        score = @x[row][seq_index] - predecessor_score(ptr)
        Step.new(:sequence_gap, nil, nil, seq_index - 1, score)
      end
    end

    def predecessor_score(ptr)
      pred_state, pred_row, pred_seq = ptr
      state_table(pred_state)[pred_row][pred_seq]
    end
  end

  # Aligns +sequence+ against +graph+.
  #
  #   result = Bio::DAGAlignment.align(graph, 'ACGT')
  #   result.score
  #   result.trace
  #
  # Options (all optional except +graph+ and +sequence+):
  # +mode+:: +:global+ (default) or +:local+
  # +match+:: score for equal characters (default 2)
  # +mismatch+:: score for different characters (default -1)
  # +gap_open+:: first-column gap cost (default -2)
  # +gap_extend+:: continuation gap cost (default -1)
  # +substitution+:: optional block receiving two characters and
  #                  returning their substitution score; when given,
  #                  +match+/+mismatch+ are ignored.
  def self.align(graph, sequence, mode: :global,
                 match: DEFAULT_MATCH, mismatch: DEFAULT_MISMATCH,
                 gap_open: DEFAULT_GAP_OPEN, gap_extend: DEFAULT_GAP_EXTEND,
                 &substitution)
    Aligner.new(graph, sequence, mode, match, mismatch,
                gap_open, gap_extend, &substitution).run
  end

end
