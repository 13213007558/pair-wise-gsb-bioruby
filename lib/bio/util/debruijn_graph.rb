#
# bio/util/debruijn_graph.rb - Bidirectional de Bruijn graph and unitig assembly
#
# Author::    BioRuby project
# License::   The Ruby License
#
# = Description
#
# Bio::DeBruijnGraph builds a bidirectional (double-strand aware) de
# Bruijn graph from short reads and compresses maximal non-branching
# paths into unitigs, without relying on an external assembler.
#
# K-mers are stored in canonical form (the lexicographically smaller
# of the k-mer and its reverse complement), so a read and its reverse
# complement contribute to exactly the same graph and never produce
# duplicate contigs. Ambiguous bases (anything other than A/C/G/T)
# break the sliding window. Isolated cycles are linearized from a
# canonically chosen start node, and every unitig is reported in a
# canonical orientation, which makes the output stable even for
# palindromic k-mers (possible when k is even) and independent of the
# order in which reads were added.
#

# = Usage
#
#   graph = Bio::DeBruijnGraph.new(31)
#   reads.each { |read| graph.add(read) }
#   graph.unitigs.each do |unitig|
#     puts "#{unitig.sequence} coverage=#{unitig.coverage}"
#   end
#
module Bio

  # Bidirectional de Bruijn graph with k-mer coverage counts.
  class DeBruijnGraph

    # A maximal non-branching path (or an isolated cycle) of the
    # graph, reported in a canonical orientation.
    class Unitig
      # Sequence of the unitig (canonical orientation).
      attr_reader :sequence
      # K-mer counts along the sequence, aligned with +sequence+.
      attr_reader :counts
      # Canonical form of the (k-1)-mer node at the left end.
      attr_reader :left_node
      # Canonical form of the (k-1)-mer node at the right end.
      attr_reader :right_node
      # Orientation of the left end node relative to its canonical
      # form: +:fwd+ or +:rev+.
      attr_reader :left_orientation
      # Orientation of the right end node relative to its canonical
      # form: +:fwd+ or +:rev+.
      attr_reader :right_orientation

      def initialize(sequence, counts, circular, k) # :nodoc:
        @sequence = sequence
        @counts = counts
        @circular = circular
        km1 = k - 1
        left = sequence[0, km1]
        right = sequence[sequence.length - km1, km1]
        @left_node = DeBruijnGraph.canonical(left)
        @left_orientation = (left == @left_node) ? :fwd : :rev
        @right_node = DeBruijnGraph.canonical(right)
        @right_orientation = (right == @right_node) ? :fwd : :rev
      end

      # Length of the unitig sequence in bases.
      def length
        @sequence.length
      end

      # Number of k-mers contained in the unitig.
      def kmer_count
        @counts.size
      end

      # Sum of the k-mer counts of the unitig.
      def total_count
        @counts.inject(0) { |sum, c| sum + c }
      end

      # Mean k-mer coverage of the unitig.
      def coverage
        total_count.to_f / @counts.size
      end

      # Whether the unitig is an isolated cycle (its last k-1 bases
      # overlap its first k-1 bases).
      def circular?
        @circular
      end
    end

    # K-mer size.
    attr_reader :k

    # Creates an empty graph for the given k-mer size.
    # _k_ must be an Integer >= 2.
    def initialize(k)
      unless k.is_a?(Integer) && k >= 2
        raise ArgumentError, "k must be an Integer greater than or equal to 2"
      end
      @k = k
      @kmer_counts = Hash.new(0)
      @built = false
    end

    # Reverse complement of a DNA string. Only A/C/G/T (both cases)
    # are complemented; other characters are left as-is.
    def self.reverse_complement(str)
      str.tr("ACGTacgt", "TGCAtgca").reverse
    end

    # Canonical form of a DNA string: the lexicographically smaller
    # of the string and its reverse complement.
    def self.canonical(str)
      rc = reverse_complement(str)
      str <= rc ? str : rc
    end

    # Adds a sequence (String or any object responding to to_s, such
    # as Bio::Sequence::NA) to the graph. Ambiguous bases break the
    # sliding k-mer window. Returns self.
    def add(sequence)
      seq = sequence.to_s.upcase
      segment_start = nil
      seq.each_char.with_index do |ch, i|
        if ch == "A" || ch == "C" || ch == "G" || ch == "T"
          segment_start = i unless segment_start
        elsif segment_start
          add_segment(seq, segment_start, i)
          segment_start = nil
        end
      end
      add_segment(seq, segment_start, seq.length) if segment_start
      @built = false
      self
    end

    # Same as DeBruijnGraph#add.
    def <<(sequence)
      add(sequence)
    end

    # Coverage count of the given k-mer. The k-mer is canonicalized
    # before lookup, so querying a k-mer or its reverse complement
    # returns the same value. Returns 0 for k-mers of the wrong
    # length or containing non-ACGT bases.
    def kmer_count(kmer)
      kmer = kmer.to_s.upcase
      return 0 unless kmer.length == @k && kmer =~ /\A[ACGT]+\z/
      @kmer_counts[self.class.canonical(kmer)]
    end

    # A Hash of canonical k-mers to coverage counts.
    def kmer_counts
      @kmer_counts.dup
    end

    # Number of canonical k-mers (edges) in the graph.
    def edge_count
      @kmer_counts.size
    end

    # Number of canonical (k-1)-mers (nodes) in the graph.
    def node_count
      build_graph unless @built
      @nodes.size
    end

    # Compresses the graph into unitigs. Returns an Array of
    # DeBruijnGraph::Unitig sorted by sequence. The result is
    # independent of the order in which reads were added.
    def unitigs
      build_graph unless @built
      visited = {}
      result = []
      sorted_edge_keys.each do |key|
        each_oriented_edge(key) do |estr|
          break if visited[key]
          source = estr[0, @k - 1]
          next if extendable?(source)
          result << walk(estr, visited)
        end
      end
      # Remaining unvisited edges form isolated cycles.
      sorted_edge_keys.each do |key|
        next if visited[key]
        result << walk_cycle(key, visited)
      end
      result.sort_by!(&:sequence)
      result
    end

    private

    def add_segment(seq, first, last)
      return if last - first < @k
      first.upto(last - @k) do |i|
        kmer = seq[i, @k]
        @kmer_counts[self.class.canonical(kmer)] += 1
      end
    end

    def sorted_edge_keys
      @kmer_counts.keys.sort
    end

    # Yields each oriented form of the canonical k-mer +key+ as a
    # String. Palindromic k-mers (only possible for even k) have a
    # single oriented form.
    def each_oriented_edge(key)
      yield key
      rc = self.class.reverse_complement(key)
      yield rc unless rc == key
    end

    # Builds the node link tables. For each canonical node (a
    # (k-1)-mer) and each of its two orientations (:fwd / :rev),
    # keeps the list of outgoing oriented edges as [key, direction]
    # pairs, where direction :fwd means the canonical k-mer itself
    # and :rc means its reverse complement.
    def build_graph
      @nodes = {}
      km1 = @k - 1
      @kmer_counts.each_key do |key|
        add_link(key[0, km1], key, :fwd)
        rc = self.class.reverse_complement(key)
        next if rc == key
        add_link(self.class.reverse_complement(key[1, km1]), key, :rc)
      end
      @built = true
    end

    def add_link(oriented_node, key, direction)
      canonical_node = self.class.canonical(oriented_node)
      node = (@nodes[canonical_node] ||= { fwd: [], rev: [] })
      node[oriented_node == canonical_node ? :fwd : :rev] << [key, direction]
    end

    # Outgoing [key, direction] links of the oriented node +u+
    # (a (k-1)-mer String in a specific orientation).
    def links_from(u)
      canonical_node = self.class.canonical(u)
      node = @nodes[canonical_node]
      return [] unless node
      node[u == canonical_node ? :fwd : :rev]
    end

    def out_degree(u)
      links_from(u).size
    end

    # An oriented node is extendable when it has exactly one incoming
    # and one outgoing edge; by the bidirectional symmetry the
    # in-degree of u equals the out-degree of its reverse complement.
    def extendable?(u)
      out_degree(u) == 1 &&
        out_degree(self.class.reverse_complement(u)) == 1
    end

    def oriented_edge_string(key, direction)
      direction == :fwd ? key : self.class.reverse_complement(key)
    end

    # Walks a maximal non-branching path starting with the oriented
    # k-mer +first_estr+, marking edges as visited.
    def walk(first_estr, visited)
      key = self.class.canonical(first_estr)
      visited[key] = true
      seq = first_estr.dup
      counts = [@kmer_counts[key]]
      target = first_estr[1, @k - 1]
      while extendable?(target)
        link = links_from(target).first
        key2 = link[0]
        break if visited[key2]
        visited[key2] = true
        estr2 = oriented_edge_string(key2, link[1])
        seq << estr2[@k - 1, 1]
        counts << @kmer_counts[key2]
        target = estr2[1, @k - 1]
      end
      make_unitig(seq, counts, false)
    end

    # Linearizes an isolated cycle containing the canonical k-mer
    # +start_key+. The cycle is reported starting from the
    # lexicographically smallest oriented node, which makes the
    # representation independent of where the traversal began.
    def walk_cycle(start_key, visited)
      km1 = @k - 1
      start_node = start_key[0, km1]
      nodes = []
      estr = start_key
      loop do
        key = self.class.canonical(estr)
        visited[key] = true
        nodes << estr[0, km1]
        target = estr[1, km1]
        break if target == start_node
        link = links_from(target).first
        estr = oriented_edge_string(link[0], link[1])
      end
      first_node = nodes.min
      seq = first_node.dup
      counts = []
      node = first_node
      begin
        link = links_from(node).first
        estr = oriented_edge_string(link[0], link[1])
        counts << @kmer_counts[link[0]]
        seq << estr[@k - 1, 1]
        node = estr[1, km1]
      end until node == first_node
      make_unitig(seq, counts, true)
    end

    # Reports every unitig in a canonical orientation (the
    # lexicographically smaller of the sequence and its reverse
    # complement), so that the output never depends on an arbitrary
    # traversal direction.
    def make_unitig(seq, counts, circular)
      rc = self.class.reverse_complement(seq)
      if rc < seq
        seq = rc
        counts = counts.reverse
      end
      Unitig.new(seq, counts, circular, @k)
    end
  end
end
