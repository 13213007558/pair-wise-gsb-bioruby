#
# = bio/tree.rb - phylogenetic tree data structure class
#
# Copyright::   Copyright (C) 2006
#               Naohisa Goto <ng@bioruby.org>
# License::     The Ruby License
#
#

require 'matrix'
require 'bio/pathway'

module Bio

  # This is the class for phylogenetic tree.
  # It stores a phylogenetic tree.
  # 
  # Internally, it is based on Bio::Pathway class.
  # However, users cannot handle Bio::Pathway object directly.
  #
  # This is alpha version. Incompatible changes may be made frequently.
  class Tree

    # Error when there are no path between specified nodes
    class NoPathError < RuntimeError; end

    # A split (bipartition) or rooted clade shared by a set of leaves.
    class Split
      include Comparable

      # Creates a new split from leaf labels.
      def initialize(labels, edge = nil)
        @labels = labels.to_a.uniq.sort.freeze
        @edge = edge
        @distance = edge && edge.respond_to?(:distance) ? edge.distance : nil
      end

      # Labels on this side of the split.
      attr_reader :labels

      # Edge represented by this split, if available.
      attr_reader :edge

      # Branch length of the represented edge.
      attr_reader :distance

      # Compares two splits by their leaf labels.
      def <=>(other)
        return nil unless other.is_a?(self.class)

        labels <=> other.labels
      end

      # Checks whether two splits contain the same leaf labels.
      def ==(other)
        other.is_a?(self.class) && labels == other.labels
      end
      alias eql? ==

      # Hash value compatible with eql?.
      def hash
        labels.hash
      end

      # Returns the labels as an array.
      def to_a
        labels.dup
      end

      # Readable representation of the split.
      def inspect
        "#<#{self.class}: #{labels.inspect}>"
      end
    end

    # Robinson–Foulds distance between two phylogenetic trees.
    class RobinsonFoulds
      # Creates a new comparison and restricts both trees to their common,
      # uniquely named leaves.
      def initialize(tree1, tree2, rooted: false, branch_length: false)
        @tree1 = tree1
        @tree2 = tree2
        @rooted = rooted
        @branch_length = branch_length

        @leaf_nodes1 = named_leaves(tree1)
        @leaf_nodes2 = named_leaves(tree2)
        leaves1 = @leaf_nodes1
        leaves2 = @leaf_nodes2
        names1 = leaves1.keys
        names2 = leaves2.keys
        @common_leaves = (names1 & names2).sort
        if @common_leaves.size < 2
          raise ArgumentError, 'trees must share at least two named leaves'
        end

        analysis1 = analyze(tree1)
        analysis2 = analyze(tree2)
        @splits1 = analysis1[:splits]
        @splits2 = analysis2[:splits]
        @first_tree_splits = (@splits1 - @splits2).sort
        @second_tree_splits = (@splits2 - @splits1).sort
        @distance = @first_tree_splits.size + @second_tree_splits.size

        calculate_branch_lengths if branch_length
      end

      # First tree in the comparison.
      attr_reader :tree1

      # Second tree in the comparison.
      attr_reader :tree2

      # Sorted names of leaves shared by both trees.
      attr_reader :common_leaves

      # Splits found in the first tree.
      attr_reader :splits1

      # Splits found in the second tree.
      attr_reader :splits2

      # Splits found only in the first tree.
      attr_reader :first_tree_splits

      # Splits found only in the second tree.
      attr_reader :second_tree_splits

      # Symmetric-difference Robinson–Foulds distance.
      attr_reader :distance

      # Absolute branch-length differences for splits shared by both trees.
      attr_reader :branch_length_differences

      # Sum of branch lengths over symmetric-difference splits.
      attr_reader :branch_length_distance

      private

      def named_leaves(tree)
        raise ArgumentError, 'rooted comparison requires tree.root' if
          @rooted && tree.root.nil?

        leaves = tree.nodes.find_all do |node|
          tree.out_degree(node) == 1
        end
        named = {}
        leaves.each do |node|
          name = tree.get_node_name(node).to_s
          next if name.empty? && node == tree.root
          if name.empty?
            raise ArgumentError, 'all leaves must have non-empty names'
          end
          if named.key?(name)
            raise ArgumentError, "duplicate leaf name: #{name}"
          end
          named[name] = node
        end
        named
      end

      def analyze(tree)
        leaf_nodes = tree.equal?(@tree1) ? @leaf_nodes1 : @leaf_nodes2
        selected = @common_leaves.map { |name| leaf_nodes.fetch(name) }
        induced = tree.subtree_with_all_paths(selected)
        adjacency = build_adjacency(induced)
        root = @rooted ? rooted_anchor(tree, adjacency) : nil
        contract_degree_two(tree, adjacency, root)
        splits = @rooted ? rooted_splits(tree, adjacency, root) :
                           unrooted_splits(tree, adjacency)
        { root: root, adjacency: adjacency, splits: splits }
      end

      def build_adjacency(tree)
        adjacency = Hash.new { |hash, node| hash[node] = [] }
        tree.edges.each do |node1, node2, edge|
          adjacency[node1] << [node2, edge]
          adjacency[node2] << [node1, edge]
        end
        adjacency
      end

      def rooted_anchor(tree, adjacency)
        root = tree.root
        raise ArgumentError, 'rooted comparison requires tree.root' unless root
        return root if adjacency.key?(root)

        target = adjacency.keys.first
        tree.path(root, target).find { |node| adjacency.key?(node) }
      end

      def contract_degree_two(tree, adjacency, root)
        pending = adjacency.keys.select do |node|
          adjacency[node].size == 2 && node != root
        end
        until pending.empty?
          node = pending.shift
          next unless adjacency.key?(node)
          neighbors = adjacency[node]
          next unless neighbors && neighbors.size == 2

          node1, edge1 = neighbors[0]
          node2, edge2 = neighbors[1]
          merged_edge = tree.get_edge_merged(edge1, edge2)
          adjacency[node1].delete_if { |other, _| other == node }
          adjacency[node2].delete_if { |other, _| other == node }
          adjacency.delete(node)
          adjacency[node1] << [node2, merged_edge]
          adjacency[node2] << [node1, merged_edge]
          [node1, node2].each do |neighbor|
            if adjacency[neighbor].size == 2 && neighbor != root
              pending << neighbor unless pending.include?(neighbor)
            end
          end
        end
      end

      def rooted_splits(tree, adjacency, root)
        splits = []
        collect_rooted_splits(tree, adjacency, root, nil, nil, splits)
        splits
      end

      def collect_rooted_splits(tree, adjacency, node, parent, edge, splits)
        children = adjacency[node].reject { |child, _| child == parent }
        name = tree.get_node_name(node).to_s
        is_leaf = adjacency[node].size == 1 && @common_leaves.include?(name)
        own_labels = is_leaf ? [name] : []
        child_labels = children.flat_map do |child, child_edge|
          collect_rooted_splits(tree, adjacency, child, node,
                                child_edge, splits)
        end
        labels = (own_labels + child_labels).sort
        if edge && labels.size >= 2 && labels.size < @common_leaves.size
          splits << Split.new(labels, edge)
        end
        labels
      end

      def unrooted_splits(tree, adjacency)
        splits = []
        visited = {}
        adjacency.each do |node1, neighbors|
          neighbors.each do |node2, edge|
            key = [node1.object_id, node2.object_id].sort
            next if visited[key]
            visited[key] = true

            side = component_labels(tree, adjacency, node1, node2)
            other = @common_leaves - side
            next unless side.size.between?(2, @common_leaves.size - 2)
            side = other if (other <=> side) < 0
            splits << Split.new(side, edge)
          end
        end
        splits
      end

      def component_labels(tree, adjacency, start, blocked)
        pending = [start]
        visited = { blocked => true, start => true }
        labels = []
        until pending.empty?
          node = pending.shift
          name = tree.get_node_name(node).to_s
          labels << name if @common_leaves.include?(name)
          adjacency[node].each do |neighbor, _|
            next if visited[neighbor]
            visited[neighbor] = true
            pending << neighbor
          end
        end
        labels.sort
      end

      def calculate_branch_lengths
        lengths1 = split_lengths(@splits1)
        lengths2 = split_lengths(@splits2)
        @branch_length_differences = {}
        (@splits1 & @splits2).each do |split|
          @branch_length_differences[split] =
            (lengths1.fetch(split, 0.0) - lengths2.fetch(split, 0.0)).abs
        end
        @branch_length_distance =
          (@splits1 | @splits2).sum do |split|
            (lengths1.fetch(split, 0.0) - lengths2.fetch(split, 0.0)).abs
          end
      end

      def split_lengths(splits)
        splits.each_with_object({}) do |split, lengths|
          lengths[split] = split.distance || 0.0
        end
      end
    end

    # Edge object of each node.
    # By default, the object doesn't contain any node information.
    class Edge

      # creates a new edge.
      def initialize(distance = nil)
        if distance.kind_of?(Numeric)
          self.distance = distance
        elsif distance
          self.distance_string = distance
        end
      end

      # evolutionary distance
      attr_reader :distance

      # evolutionary distance represented as a string
      attr_reader :distance_string

      # set evolutionary distance value
      def distance=(num)
        @distance = num
        @distance_string = (num ? num.to_s : num)
      end

      # set evolutionary distance value from a string
      def distance_string=(str)
        if str.to_s.strip.empty?
          @distance = nil
          @distance_string = str
        else
          @distance = str.to_f
          @distance_string = str
        end
      end

      # visualization of this object
      def inspect
        "<Edge distance=#{@distance.inspect}>"
      end

      # string representation of this object
      def to_s
        @distance_string.to_s
      end

      #---
      # methods for NHX (New Hampshire eXtended) and/or PhyloXML
      #+++

      # log likelihood value (:L in NHX)
      attr_accessor :log_likelihood

      # width of the edge
      # (<branch width="w"> of PhyloXML, or :W="w" in NHX)
      attr_accessor :width

      # Other NHX parameters. Returns a Hash.
      # Note that :L and :W
      # are not stored here but stored in the proper attributes in this class.
      # However, if you force to set these parameters in this hash,
      # the parameters in this hash are preferred when generating NHX.
      # In addition, If the same parameters are defined at Node object,
      # the parameters in the node are preferred.
      def nhx_parameters
        @nhx_parameters ||= {}
        @nhx_parameters
      end

    end #class Edge

    # Gets distance value from the given edge.
    # Returns float or any other numeric value or nil.
    def get_edge_distance(edge)
      begin
        dist = edge.distance
      rescue NoMethodError
        dist = edge
      end
      dist
    end

    # Gets distance string from the given edge.
    # Returns a string or nil.
    def get_edge_distance_string(edge)
      begin
        dist = edge.distance_string
      rescue NoMethodError
        dist = (edge ? edge.to_s : nil)
      end
      dist
    end

    # Returns edge1 + edge2
    def get_edge_merged(edge1, edge2)
      dist1 = get_edge_distance(edge1)
      dist2 = get_edge_distance(edge2)
      if dist1 and dist2 then
        Edge.new(dist1 + dist2)
      elsif dist1 then
        Edge.new(dist1)
      elsif dist2 then
        Edge.new(dist2)
      else
        Edge.new
      end
    end

    # Node object.
    class Node

      # Creates a new node.
      def initialize(name = nil)
        @name = name if name
      end

      # name of the node
      attr_accessor :name

      # bootstrap value
      attr_reader :bootstrap

      # bootstrap value as a string
      attr_reader :bootstrap_string

      # sets a bootstrap value
      def bootstrap=(num)
        @bootstrap_string = (num ? num.to_s : num)
        @bootstrap = num
      end

      # sets a bootstrap value from a string
      def bootstrap_string=(str)
        if str.to_s.strip.empty?
          @bootstrap = nil
          @bootstrap_string = str
        else
          i = str.to_i
          f = str.to_f
          @bootstrap = (i == f ? i : f)
          @bootstrap_string = str
        end
      end

      # visualization of this object
      def inspect
        if @name and !@name.empty? then
          str = "(Node:#{@name.inspect}"
        else
          str = sprintf('(Node:%x', (self.__id__ << 1) & 0xffffffff)
        end
        if defined?(@bootstrap) and @bootstrap then
          str += " bootstrap=#{@bootstrap.inspect}"
        end
        str += ")"
        str
      end

      # string representation of this object
      def to_s
        @name.to_s
      end

      # the order of the node
      # (lower value, high priority)
      attr_accessor :order_number

      #---
      # methods for NHX (New Hampshire eXtended) and/or PhyloXML
      #+++

      # Phylogenetic events.
      # Returns an Array of one (or more?) of the following symbols
      #   :gene_duplication
      #   :speciation
      def events
        @events ||= []
        @events
      end

      # EC number (EC_number in PhyloXML, or :E in NHX)
      attr_accessor :ec_number

      # scientific name (scientific_name in PhyloXML, or :S in NHX)
      attr_accessor :scientific_name

      # taxonomy identifier (taxonomy_identifier in PhyloXML, or :T in NHX)
      attr_accessor :taxonomy_id

      # Other NHX parameters. Returns a Hash.
      # Note that :D, :E, :S, and :T
      # are not stored here but stored in the proper attributes in this class.
      # However, if you force to set these parameters in this hash,
      # the parameters in this hash are preferred when generating NHX.
      def nhx_parameters
        @nhx_parameters ||= {}
        @nhx_parameters
      end

    end #class Node

    # Gets node name
    def get_node_name(node)
      begin
        node.name
      rescue NoMethodError
        node.to_s
      end
    end

    def get_node_bootstrap(node)
      begin
        node.bootstrap
      rescue NoMethodError
        nil
      end
    end

    def get_node_bootstrap_string(node)
      begin
        node.bootstrap_string
      rescue NoMethodError
        nil
      end
    end

    # Creates a new phylogenetic tree.
    # When no arguments are given, it creates a new empty tree.
    # When a Tree object is given, it copies the tree.
    # Note that the  new tree shares Node and Edge objects
    # with the given tree.
    def initialize(tree = nil)
      # creates an undirected adjacency list graph
      @pathway = Bio::Pathway.new([], true)
      @root = nil
      @options = {}
      _init_cache
      self.concat(tree) if tree
    end

    # (private) clear internal cache
    def _init_cache
      @cache_parent = {}
    end
    private :_init_cache

    # (private) clear internal cache
    def _clear_cache
      @cache_parent.clear
    end
    private :_clear_cache

    # root node of this tree
    # (even if unrooted tree, it is used by some methods)
    attr_accessor :root

    # tree options; mainly used for tree output
    attr_accessor :options

    # Clears all nodes and edges.
    # Returns self.
    # Note that options and root are also cleared.
    def clear
      initialize
      self
    end

    # Returns all nodes as an array.
    def nodes
      @pathway.graph.keys
    end

    # Number of nodes.
    def number_of_nodes
      @pathway.nodes
    end

    # Iterates over each node of this tree.
    def each_node(&x) #:yields: node
      @pathway.graph.each_key(&x)
      self
    end

    # Iterates over each edges of this tree.
    def each_edge #:yields: source, target, edge
      @pathway.relations.each do |rel|
        yield rel.node[0], rel.node[1], rel.relation
      end
      self
    end

    # Returns all edges an array of [ node0, node1, edge ]
    def edges
      @pathway.relations.collect do |rel|
        [ rel.node[0], rel.node[1], rel.relation ]
      end
    end

    # Returns number of edges in the tree.
    def number_of_edges
      @pathway.relations.size
    end

    # Returns an array of adjacent nodes of the given node.
    def adjacent_nodes(node)
      h = @pathway.graph[node]
      h ? h.keys : []
    end

    # Returns all connected edges with adjacent nodes.
    # Returns an array of the array [ source, target, edge ].
    #
    # The reason why the method name is "out_edges" is that
    # it comes from the Boost Graph Library.
    def out_edges(source)
      h = @pathway.graph[source]
      if h
        h.collect { |key, val| [ source, key, val ] }
      else
        []
      end
    end

    # Iterates over each connected edges of the given node.
    # Returns self.
    #
    # The reason why the method name is "each_out_edge" is that
    # it comes from the Boost Graph Library.
    def each_out_edge(source) #:yields: source, target, edge
      h = @pathway.graph[source]
      h.each { |key, val| yield source, key, val } if h
      self
    end

    # Returns number of edges in the given node.
    #
    # The reason why the method name is "out_degree" is that
    # it comes from the Boost Graph Library.
    def out_degree(source)
      h = @pathway.graph[source]
      h ? h.size : 0
    end

    # Returns an edge from source to target.
    # If source and target are not adjacent nodes, returns nil.
    def get_edge(source, target)
      h = @pathway.graph[source]
      h ? h[target] : nil
    end

    # Adds a new edge to the tree.
    # Returns the newly added edge.
    # If the edge already exists, it is overwritten with new one.
    def add_edge(source, target, edge = Edge.new)
      _clear_cache
      @pathway.append(Bio::Relation.new(source, target, edge))
      edge
    end

    # Finds a node in the tree by given name and returns the node.
    # If the node does not found, returns nil.
    # If multiple nodes with the same name exist,
    # the result would be one of those (unspecified).
    def get_node_by_name(str)
      self.each_node do |node|
        if get_node_name(node) == str
          return node
        end
      end
      nil
    end

    # Adds a node to the tree.
    # Returns self.
    # If the node already exists, it does nothing.
    def add_node(node)
      _clear_cache
      @pathway.graph[node] ||= {}
      self
    end

    # If the node exists, returns true.
    # Otherwise, returns false.
    def include?(node)
      @pathway.graph[node] ? true : false
    end

    # Removes all edges connected with the node.
    # Returns self.
    # If the node does not exist, raises IndexError.
    def clear_node(node)
      unless self.include?(node)
        raise IndexError, 'the node does not exist'
      end
      _clear_cache
      @pathway.relations.delete_if do |rel|
        rel.node.include?(node)
      end
      @pathway.graph[node].each_key do |k|
        @pathway.graph[k].delete(node)
      end
      @pathway.graph[node].clear
      self
    end

    # Removes the given node from the tree.
    # All edges connected with the node are also removed.
    # Returns self.
    # If the node does not exist, raises IndexError.
    def remove_node(node)
      #_clear_cache #done in clear_node(node)
      self.clear_node(node)
      @pathway.graph.delete(node)
      self
    end

    # Removes each node if the block returns not nil.
    # All edges connected with the removed nodes are also removed.
    # Returns self.
    def remove_node_if
      #_clear_cache #done in clear_node(node)
      all = self.nodes
      all.each do |node|
        if yield node then
          self.clear_node(node)
          @pathway.graph.delete(node)
        end
      end
      self
    end

    # Removes an edge between source and target.
    # Returns self.
    # If the edge does not exist, raises IndexError.
    #---
    # If two or more edges exists between source and target,
    # all of them are removed.
    #+++
    def remove_edge(source, target)
      unless self.get_edge(source, target) then
        raise IndexError, 'edge not found'
      end
      _clear_cache
      fwd = [ source, target ]
      rev = [ target, source ]
      @pathway.relations.delete_if do |rel|
        rel.node == fwd or rel.node == rev
      end
      h = @pathway.graph[source]
      h.delete(target) if h
      h = @pathway.graph[target]
      h.delete(source) if h
      self
    end

    # Removes each edge if the block returns not nil.
    # Returns self.
    def remove_edge_if #:yields: source, target, edge
      _clear_cache
      removed_rel = []
      @pathway.relations.delete_if do |rel|
        if yield rel.node[0], rel.node[1], rel.edge then
          removed_rel << rel
          true
        end
      end
      removed_rel.each do |rel|
        source = rel.node[0]
        target = rel.node[1]
        h = @pathway.graph[source]
        h.delete(target) if h
        h = @pathway.graph[target]
        h.delete(source) if h
      end
      self
    end

    # Replaces each node by each block's return value.
    # Returns self.
    def collect_node! #:yields: node
      _clear_cache
      tr = {}
      self.each_node do |node|
        tr[node] = yield node
      end
      # replaces nodes in @pathway.relations
      @pathway.relations.each do |rel|
        rel.node.collect! { |node| tr[node] }
      end
      # re-generates @pathway from relations
      @pathway.to_list
      # adds orphan nodes
      tr.each_value do |newnode|
        @pathway.graph[newnode] ||= {}
      end
      self
    end

    # Replaces each edge by each block's return value.
    # Returns self.
    def collect_edge! #:yields: source, target, edge
      _clear_cache
      @pathway.relations.each do |rel|
        newedge = yield rel.node[0], rel.node[1], rel.relation
        rel.edge = newedge
        @pathway.append(rel, false)
      end
      self
    end

    # Gets the sub-tree consisted of given nodes.
    # _nodes_ must be an array of nodes.
    # Nodes that do not exist in the original tree are ignored.
    # Returns a Tree object.
    # Note that the sub-tree shares Node and Edge objects
    # with the original tree.
    def subtree(nodes)
      nodes = nodes.find_all do |x|
        @pathway.graph[x]
      end
      return self.class.new if nodes.empty?
      # creates subtree
      new_tree = self.class.new
      nodes.each do |x|
        new_tree.add_node(x)
      end
      self.each_edge do |node1, node2, edge|
        if new_tree.include?(node1) and new_tree.include?(node2) then
          new_tree.add_edge(node1, node2, edge)
        end
      end
      return new_tree
    end

    # Gets the sub-tree consisted of given nodes and
    # all internal nodes connected between given nodes.
    # _nodes_ must be an array of nodes.
    # Nodes that do not exist in the original tree are ignored.
    # Returns a Tree object.
    # The result is unspecified for cyclic trees.
    # Note that the sub-tree shares Node and Edge objects
    # with the original tree.
    def subtree_with_all_paths(nodes)
      hash = {}
      nodes.each { |x| hash[x] = true }
      nodes.each_index do |i|
        node1 = nodes[i]
        (0...i).each do |j|
          node2 = nodes[j]
          unless node1 == node2 then
            begin
              path = self.path(node1, node2)
            rescue IndexError, NoPathError
              path = []
            end
            path.each { |x| hash[x] = true }
          end
        end
      end
      self.subtree(hash.keys)
    end

    # Concatenates the other tree.
    # If the same edge exists, the edge in _other_ is used.
    # Returns self.
    # The result is unspecified if _other_ isn't a Tree object.
    # Note that the Node and Edge objects in the _other_ tree are
    # shared in the concatinated tree.
    def concat(other)
      #raise TypeError unless other.kind_of?(self.class)
      _clear_cache
      other.each_node do |node|
        self.add_node(node)
      end
      other.each_edge do |node1, node2, edge|
        self.add_edge(node1, node2, edge)
      end
      self
    end

    # Gets path from node1 to node2.
    # Returns an array of nodes, including node1 and node2.
    # If node1 and/or node2 do not exist, IndexError is raised.
    # If node1 and node2 are not connected, NoPathError is raised.
    # The result is unspecified for cyclic trees.
    def path(node1, node2)
      raise IndexError, 'node1 not found' unless @pathway.graph[node1]
      raise IndexError, 'node2 not found' unless @pathway.graph[node2]
      return [ node1 ] if node1 == node2
      return [ node1, node2 ] if @pathway.graph[node1][node2]
      _, path = @pathway.bfs_shortest_path(node1, node2)
      unless path[0] == node1 and path[-1] == node2 then
        raise NoPathError, 'node1 and node2 are not connected'
      end
      path
    end

    # Iterates over each edge from node1 to node2.
    # The result is unspecified for cyclic trees.
    def each_edge_in_path(node1, node2)
      path = self.path(node1, node2)
      source = path.shift
      path.each do |target|
        edge = self.get_edge(source, target)
        yield source, target, edge
        source = target
      end
      self
    end

    # Returns distance between node1 and node2.
    # It would raise error if the edges didn't contain distance values.
    # The result is unspecified for cyclic trees.
    def distance(node1, node2)
      distance = 0
      self.each_edge_in_path(node1, node2) do |source, target, edge|
        distance += get_edge_distance(edge)
      end
      distance
    end

    # (private) get parent only by using cache
    def _get_cached_parent(node, root)
      @cache_parent[root] ||= Hash.new
      cache = @cache_parent[root]
      if node == root then
        unless cache.has_key?(root) then
          self.adjacent_nodes(root).each do |n|
            cache[n] ||= root if n != root
          end
          cache[root] = nil
        end
        parent = nil
      else
        unless parent = cache[node] then
          parent = self.adjacent_nodes(node).find { |n|
            (m = cache[n]) && (m != node)
          }
          _cache_parent(node, parent, root) if parent
        end
        parent
      end
    end
    private :_get_cached_parent

    # (private) set parent cache
    def _cache_parent(node, parent, root)
      return unless parent
      cache = @cache_parent[root]
      cache[node] = parent
      self.adjacent_nodes(node).each do |n|
        cache[n] ||= node if n != parent
      end
    end
    private :_cache_parent

    # Gets the parent node of the _node_.
    # If _root_ isn't specified or _root_ is <code>nil</code>, @root is used.
    # Returns an <code>Node</code> object or nil.
    # The result is unspecified for cyclic trees.
    def parent(node, root = nil)
      root ||= @root
      raise IndexError, 'can not get parent for unrooted tree' unless root
      unless ret = _get_cached_parent(node, root) then
        ret = self.path(root, node)[-2]
        _cache_parent(node, ret, root)
      end
      ret
    end

    # Gets the adjacent children nodes of the _node_.
    # If _root_ isn't specified or _root_ is <code>nil</code>, @root is used.
    # Returns an array of <code>Node</code>s.
    # The result is unspecified for cyclic trees.
    def children(node, root = nil)
      root ||= @root
      c = self.adjacent_nodes(node)
      c.delete(self.parent(node, root))
      c
    end

    # Gets all descendent nodes of the _node_.
    # If _root_ isn't specified or _root_ is <code>nil</code>, @root is used.
    # Returns an array of <code>Node</code>s.
    # The result is unspecified for cyclic trees.
    def descendents(node, root = nil)
      root ||= @root
      distance, route = @pathway.breadth_first_search(root)
      d = distance[node]
      result = []
      distance.each do |key, val|
        if val > d then
          x = key
          while x = route[x]
            if x == node then
              result << key
              break
            end
            break if distance[x] <= d
          end
        end
      end
      result
    end

    # If _node_ is nil, returns an array of 
    # all leaves (nodes connected with one edge).
    # Otherwise, gets all descendent leaf nodes of the _node_.
    # If _root_ isn't specified or _root_ is <code>nil</code>, @root is used.
    # Returns an array of <code>Node</code>s.
    # The result is unspecified for cyclic trees.
    def leaves(node = nil, root = nil)
      unless node then
        nodes = []
        self.each_node do |x|
          nodes << x if self.out_degree(x) == 1
        end
        return nodes
      else
        root ||= @root
        self.descendents(node, root).find_all do |x|
          self.adjacent_nodes(x).size == 1
        end
      end
    end

    # Gets all ancestral nodes of the _node_.
    # If _root_ isn't specified or _root_ is <code>nil</code>, @root is used.
    # Returns an array of <code>Node</code>s.
    # The result is unspecified for cyclic trees.
    def ancestors(node, root = nil)
      root ||= @root
      (self.path(root, node) - [ node ]).reverse
    end

    # Gets the lowest common ancestor of the two nodes.
    # If _root_ isn't specified or _root_ is <code>nil</code>, @root is used.
    # Returns a <code>Node</code> object or nil.
    # The result is unspecified for cyclic trees.
    def lowest_common_ancestor(node1, node2, root = nil)
      root ||= @root
      _, route = @pathway.breadth_first_search(root)
      x = node1; r1 = []
      begin; r1 << x; end while x = route[x]
      x = node2; r2 = []
      begin; r2 << x; end while x = route[x]
      return (r1 & r2).first
    end

    # Returns total distance of all edges.
    # It would raise error if some edges didn't contain distance values.
    def total_distance
      distance = 0
      self.each_edge do |source, target, edge|
        distance += get_edge_distance(edge)
      end
      distance
    end

    # Calculates distance matrix of given nodes.
    # If _nodes_ is nil, or is ommited, it acts the same as
    # tree.distance_matrix(tree.leaves).
    # Returns a matrix object.
    # The result is unspecified for cyclic trees.
    # Note 1: The diagonal values of the matrix are 0.
    # Note 2: If the distance cannot be calculated, nil will be set.
    def distance_matrix(nodes = nil)
      nodes ||= self.leaves
      matrix = []
      nodes.each_index do |i|
        row = []
        nodes.each_index do |j|
          if i == j then
            distance = 0
          elsif r = matrix[j] and val = r[i] then
            distance = val
          else
            distance = (self.distance(nodes[i], nodes[j]) rescue nil)
          end
          row << distance
        end
        matrix << row
      end
      Matrix.rows(matrix, false)
    end

    # Shows the adjacency matrix representation of the tree.
    # It shows matrix only for given nodes.
    # If _nodes_ is nil or is ommitted,
    # it acts the same as tree.adjacency_matrix(tree.nodes).
    # If a block is given, for each edge,
    # it yields _source_, _target_, and _edge_, and
    # uses the returned value of the block.
    # Without blocks, it uses edge.
    # Returns a matrix object.
    def adjacency_matrix(nodes = nil,
                         default_value = nil,
                         diagonal_value = nil) #:yields: source, target, edge
      nodes ||= self.nodes
      size = nodes.size
      hash = {}
      nodes.each_with_index { |x, i| hash[x] = i }
      # prepares an matrix
      matrix = Array.new(size, nil)
      matrix.collect! { |x| Array.new(size, default_value) }
      (0...size).each { |i| matrix[i][i] = diagonal_value }
      # fills the matrix from each edge
      self.each_edge do |source, target, edge|
        i_source = hash[source]
        i_target = hash[target]
        if i_source and i_target then
          val = block_given? ? (yield source, target, edge) : edge
          matrix[i_source][i_target] = val
          matrix[i_target][i_source] = val
        end
      end
      Matrix.rows(matrix, false)
    end

    # Removes all nodes that are not branches nor leaves.
    # That is, removes nodes connected with exactly two edges.
    # For each removed node, two adjacent edges are merged and
    # a new edge are created.
    # Returns removed nodes.
    # Note that orphan nodes are still kept unchanged.
    def remove_nonsense_nodes
      _clear_cache
      hash = {}
      self.each_node do |node|
        hash[node] = true if @pathway.graph[node].size == 2
      end
      hash.each_key do |node|
        adjs = @pathway.graph[node].keys
        edges = @pathway.graph[node].values
        new_edge = get_edge_merged(edges[0], edges[1])
        @pathway.graph[adjs[0]].delete(node)
        @pathway.graph[adjs[1]].delete(node)
        @pathway.graph.delete(node)
        @pathway.append(Bio::Relation.new(adjs[0], adjs[1], new_edge))
      end
      #@pathway.to_relations
      @pathway.relations.reject! do |rel|
        hash[rel.node[0]] or hash[rel.node[1]]
      end
      return hash.keys
    end

    # Insert a new node between adjacent nodes node1 and node2.
    # The old edge between node1 and node2 are changed to the edge
    # between new_node and node2.
    # The edge between node1 and new_node is newly created.
    #
    # If new_distance is specified, the distance between
    # node1 and new_node is set to new_distance, and
    # distance between new_node and node2 is set to
    # <code>tree.get_edge(node1, node2).distance - new_distance</code>.
    #
    # Returns self.
    # If node1 and node2 are not adjacent, raises IndexError.
    #
    # If new_node already exists in the tree, the tree would become
    # circular. In addition, if the edge between new_node and
    # node1 (or node2) already exists, it will be erased.
    def insert_node(node1, node2, new_node, new_distance = nil)
      unless edge = self.get_edge(node1, node2) then
        raise IndexError, 'nodes not found or two nodes are not adjacent'
      end
      _clear_cache
      new_edge = Edge.new(new_distance)
      self.remove_edge(node1, node2)
      self.add_edge(node1, new_node, new_edge)
      if new_distance and old_distance = get_edge_distance(edge) then
        old_distance -= new_distance
        begin
          edge.distance = old_distance
        rescue NoMethodError
          edge = old_distance
        end
      end
      self.add_edge(new_node, node2, edge)
      self
    end

    # Compares this tree with another tree by Robinson–Foulds distance.
    # Trees are first restricted to their uniquely named common leaves.
    # Set rooted: true for rooted (clade) comparison; the default is
    # unrooted (bipartition) comparison.
    # Set branch_length: true to also compare branch lengths.
    def robinson_foulds(other, rooted: false, branch_length: false)
      RobinsonFoulds.new(self, other,
                         rooted: rooted,
                         branch_length: branch_length)
    end

    # Returns the Robinson–Foulds distance from another tree.
    def robinson_foulds_distance(other, rooted = false, branch_length = false)
      comparison = robinson_foulds(other,
                                   rooted: rooted,
                                   branch_length: branch_length)
      if branch_length
        comparison.branch_length_distance
      else
        comparison.distance
      end
    end
  end #class Tree
end #module Bio

#---
# temporary added
#+++
require 'bio/tree/output'
