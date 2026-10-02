#
# = bio/tree/robinson_foulds.rb - Robinson-Foulds distances
#

module Bio
  class Tree
    # A bipartition (unrooted trees) or clade (rooted trees).
    class Split
      attr_reader :left, :right, :branch_length
      attr_writer :branch_length

      def initialize(left, right = nil, branch_length = nil)
        @left = left.dup.freeze
        @right = right ? right.dup.freeze : nil
        @branch_length = branch_length
      end

      def rooted?
        @right.nil?
      end

      def key
        rooted? ? @left : [ @left, @right ]
      end

      def names
        rooted? ? @left : [ @left, @right ]
      end

      def clade
        @left
      end

      def ==(other)
        return false unless other.is_a?(self.class)

        key == other.key
      end
      alias eql? ==

      def hash
        key.hash
      end

      def inspect
        "#<Bio::Tree::Split #{names.inspect}>"
      end
    end

    # Result of a Robinson-Foulds comparison.
    class RobinsonFouldsResult
      attr_reader :distance, :common_leaves,
                  :tree1_only_splits, :tree2_only_splits,
                  :branch_length_distance

      def initialize(distance, common_leaves, tree1_only_splits,
                     tree2_only_splits, branch_length_distance)
        @distance = distance
        @common_leaves = common_leaves
        @tree1_only_splits = tree1_only_splits
        @tree2_only_splits = tree2_only_splits
        @branch_length_distance = branch_length_distance
      end

      def only_in_tree1
        @tree1_only_splits
      end
      alias splits_only_in_tree1 only_in_tree1

      def only_in_tree2
        @tree2_only_splits
      end
      alias splits_only_in_tree2 only_in_tree2

      def to_i
        @distance
      end
      alias to_int to_i

      def to_f
        @distance.to_f
      end

      def ==(other)
        if other.is_a?(self.class)
          @distance == other.distance &&
            @common_leaves == other.common_leaves &&
            @tree1_only_splits == other.only_in_tree1 &&
            @tree2_only_splits == other.only_in_tree2
        else
          @distance == other
        end
      end
      alias eql? ==

      def hash
        @distance.hash
      end

      def inspect
        "#<Bio::Tree::RobinsonFouldsResult distance=#{@distance} " \
          "common_leaves=#{@common_leaves.inspect}>"
      end
    end

    # Calculates rooted and unrooted Robinson-Foulds distances after
    # restricting both trees to their uniquely named common leaves.
    class RobinsonFoulds
      def initialize(tree1, tree2, rooted: false, branch_lengths: false)
        @tree1 = tree1
        @tree2 = tree2
        @rooted = rooted
        @branch_lengths = branch_lengths
      end

      def self.normalize_options(args, options)
        options = options.dup
        if args.first.is_a?(Hash)
          options.merge!(args.shift)
        end
        rooted = args[0].nil? ? options.fetch(:rooted, false) : args[0]
        branch_lengths =
          if args[1].nil?
            options.fetch(:branch_lengths, false) ||
              options.fetch(:branch_length, false) ||
              options.fetch(:weighted, false)
          else
            args[1]
          end
        [ rooted, branch_lengths ]
      end

      def compare
        leaves1 = named_leaves(@tree1)
        leaves2 = named_leaves(@tree2)
        common = (leaves1.keys & leaves2.keys).sort
        validate_tree(@tree1)
        validate_tree(@tree2)
        if common.empty?
          return RobinsonFouldsResult.new(0, [], [], [],
                                          @branch_lengths ? 0.0 : nil)
        end

        pruned1 = prepared_tree(@tree1, common.map { |name| leaves1[name] })
        pruned2 = prepared_tree(@tree2, common.map { |name| leaves2[name] })
        @pruned1 = pruned1
        @pruned2 = pruned2

        splits1 = splits(pruned1, common)
        splits2 = splits(pruned2, common)

        only1 = splits1.reject { |key, _split| splits2.key?(key) }
        only2 = splits2.reject { |key, _split| splits1.key?(key) }

        distance = only1.size + only2.size
        branch_length_distance =
          if @branch_lengths
            branch_length_difference(splits1, splits2, only1, only2)
          end

        RobinsonFouldsResult.new(
          distance,
          common,
          only1.values.map { |split| split.names }.sort,
          only2.values.map { |split| split.names }.sort,
          branch_length_distance
        )
      end

      private

      def named_leaves(tree)
        leaves = {}
        tree.nodes.each do |node|
          next if tree.out_degree(node) > 1

          name = tree.get_node_name(node).to_s
          next if name.empty?
          if leaves.key?(name)
            raise ArgumentError, "duplicate leaf name: #{name}"
          end

          leaves[name] = node
        end
        leaves
      end

      def validate_tree(tree)
        connected = {}
        start = tree.nodes.first
        if start
          connected[start] = true
          queue = [ start ]
          until queue.empty?
            node = queue.shift
            tree.adjacent_nodes(node).each do |neighbor|
              next if connected.key?(neighbor)

              connected[neighbor] = true
              queue << neighbor
            end
          end
        end
        unless connected.size == tree.number_of_nodes
          raise ArgumentError, 'Robinson-Foulds requires a connected tree'
        end
      end

      def prepared_tree(tree, leaf_nodes)
        pruned = induced_tree(tree, leaf_nodes)
        if @rooted
          root = tree.root
          if root.nil?
            raise ArgumentError, 'rooted comparison requires tree.root'
          end

          unless pruned.include?(root)
            root = tree.path(root, leaf_nodes.first).find do |node|
              pruned.include?(node)
            end
          end
          suppress_degree_two(pruned, root)
          pruned.root = root
        else
          suppress_degree_two(pruned, nil)
          pruned.root = pruned.nodes.first
        end
        pruned
      end

      def induced_tree(tree, leaf_nodes)
        result = Tree.new
        return result if leaf_nodes.empty?

        kept = {}
        leaf_nodes.each { |node| kept[node] = true }

        if leaf_nodes.size > 1
          leaf_nodes.combination(2) do |node1, node2|
            tree.path(node1, node2).each { |node| kept[node] = true }
          end
        end

        result.concat(tree.subtree(kept.keys))
      end

      def suppress_degree_two(tree, preserved_node)
        loop do
          node = tree.nodes.find do |candidate|
            candidate != preserved_node && tree.out_degree(candidate) == 2
          end
          break unless node

          neighbors = tree.adjacent_nodes(node)
          merged_edge = tree.get_edge_merged(
            tree.get_edge(node, neighbors[0]),
            tree.get_edge(node, neighbors[1])
          )
          suppress_node(tree, node, neighbors[0], neighbors[1],
                         merged_edge)
        end
      end

      def suppress_node(tree, node, neighbor1, neighbor2, merged_edge)
        pathway = tree.instance_variable_get(:@pathway)
        pathway.relations.reject! do |relation|
          relation.node.include?(node)
        end
        pathway.graph[neighbor1].delete(node)
        pathway.graph[neighbor2].delete(node)
        pathway.graph.delete(node)
        pathway.graph[neighbor1][neighbor2] = merged_edge
        pathway.graph[neighbor2][neighbor1] = merged_edge
        pathway.relations << Bio::Relation.new(neighbor1, neighbor2,
                                               merged_edge)
      end

      def splits(tree, common_leaves)
        result = {}
        return result if common_leaves.size < 3

        common_set = {}
        common_leaves.each { |name| common_set[name] = true }
        leaf_nodes = {}
        tree.nodes.each do |node|
          next if tree.out_degree(node) > 1

          name = tree.get_node_name(node).to_s
          leaf_nodes[node] = name if common_set.key?(name)
        end

        clusters = descendent_leaf_clusters(tree, leaf_nodes)
        clusters.each do |node, parent, cluster|
          next unless include_cluster?(cluster, common_leaves.size)

          edge = tree.get_edge(parent, node)
          if @rooted
            split = Split.new(cluster)
          else
            complement = common_leaves - cluster
            left, right = [ cluster, complement ].sort
            split = Split.new(left, right)
          end
          split.branch_length = split_branch_length(edge)
          result[split.key] = split
        end
        result
      end

      def include_cluster?(cluster, leaf_count)
        if @rooted
          cluster.size >= 2 && cluster.size < leaf_count
        else
          cluster.size >= 2 && cluster.size <= leaf_count - 2
        end
      end

      def split_branch_length(edge)
        return nil unless @branch_lengths

        length = @pruned1&.get_edge_distance(edge)
        length ||= @pruned2&.get_edge_distance(edge)
        length || 0
      end

      def branch_length_difference(splits1, splits2, only1, only2)
        total = 0.0
        only1.each_value { |split| total += split.branch_length }
        only2.each_value { |split| total += split.branch_length }
        splits1.each do |key, split|
          next unless splits2.key?(key)

          total += (split.branch_length - splits2[key].branch_length).abs
        end
        total
      end

      def descendent_leaf_clusters(tree, leaf_nodes)
        clusters = {}
        stack = [ [ tree.root, nil, false ] ]

        until stack.empty?
          node, parent, processed = stack.pop
          if processed
            names = []
            tree.adjacent_nodes(node).each do |child|
              next if child == parent

              if leaf_nodes.key?(child)
                names << leaf_nodes[child]
              else
                names.concat(clusters[[ child, node ]])
              end
            end
            clusters[[ node, parent ]] = names.sort
            next
          end

          stack.push([ node, parent, true ])
          tree.adjacent_nodes(node).each do |child|
            next if child == parent

            stack.push([ child, node, false ])
          end
        end

        clusters.each_with_object([]) do |(key, cluster), result|
          node, parent = key
          result << [ node, parent, cluster ] unless parent.nil?
        end
      end
    end

    # Compares this tree with another tree.
    # Returns a RobinsonFouldsResult.
    def robinson_foulds(other, *args, **options)
      rooted, branch_lengths =
        RobinsonFoulds.normalize_options(args, options)
      Tree.robinson_foulds(self, other,
                           rooted: rooted,
                           branch_lengths: branch_lengths)
    end
    alias rf_distance robinson_foulds
    alias robinson_foulds_distance robinson_foulds

    class << self
      # Compares two trees on their uniquely named common leaves.
      def robinson_foulds(tree1, tree2, *args, **options)
        rooted, branch_lengths =
          RobinsonFoulds.normalize_options(args, options)
        RobinsonFoulds.new(tree1, tree2,
                           rooted: rooted,
                           branch_lengths: branch_lengths).compare
      end
      alias rf_distance robinson_foulds
      alias robinson_foulds_distance robinson_foulds
    end
  end
end
