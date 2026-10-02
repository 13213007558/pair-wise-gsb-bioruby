#
# = bio/phylo/likelihood.rb - phylogenetic likelihood of an alignment on a tree
#
# License::     The Ruby License
#

require 'bio/tree'

module Bio

  # == Description
  #
  # Bio::Phylo collects phylogenetic inference utilities that
  # operate on Bio::Tree objects.
  #
  module Phylo

    # == Description
    #
    # The Jukes-Cantor 1969 nucleotide substitution model.
    #
    # All four bases have equilibrium frequency 1/4 and all
    # substitutions occur at the same rate.  Branch lengths are
    # measured in expected numbers of substitutions per site, so the
    # transition probabilities for a branch of length +t+ are
    #
    #   P_ii(t) = 1/4 + 3/4 * exp(-4t/3)   (same state)
    #   P_ij(t) = 1/4 - 1/4 * exp(-4t/3)   (different state)
    #
    # For +t+ == 0 this is exactly the identity matrix: a zero-length
    # branch deterministically transmits the ancestral state.
    #
    class JC69

      # Nucleotide states, in the fixed order used by all vectors
      # and matrices of this class.
      STATES = %w[A C G T].freeze

      # Equilibrium state frequencies (uniform for JC69).
      # Returns an Array of 4 Floats in STATES order.
      def state_frequencies
        [0.25, 0.25, 0.25, 0.25]
      end

      # Transition probability matrix P(t) for a branch of length +t+
      # (expected number of substitutions per site; already scaled by
      # the rate category when rate mixtures are used).
      # Returns a 4x4 nested Array in STATES order:
      # entry [i][j] is the probability of state j given state i.
      def transition_probabilities(t)
        t = Float(t)
        if t < 0 then
          raise ArgumentError, "negative branch length: #{t.inspect}"
        end
        e = Math.exp(-4.0 * t / 3.0)
        diag = 0.25 + 0.75 * e
        off  = 0.25 - 0.25 * e
        [0, 1, 2, 3].map do |i|
          [0, 1, 2, 3].map { |j| i == j ? diag : off }
        end
      end
    end #class JC69

    # == Description
    #
    # IUPAC degenerate nucleotide symbols, interpreted as the set of
    # concrete states a site may take.  State indices follow
    # Bio::Phylo::JC69::STATES order (A, C, G, T).
    #
    module IUPAC

      # All four state indices (A, C, G, T).
      ALL_STATES = [0, 1, 2, 3].freeze

      # Mapping of IUPAC symbols (upper case) to allowed state
      # indices.  'U' is treated as 'T'.
      #
      # Gap policy: the gap characters '-' and '.' (and the unknown
      # character '?') are treated as completely missing data, i.e.
      # all four states are allowed, exactly like 'N'.  Such
      # observations contribute no information to the likelihood.
      TABLE = {
        'A' => [0], 'C' => [1], 'G' => [2], 'T' => [3], 'U' => [3],
        'R' => [0, 2], 'Y' => [1, 3], 'S' => [1, 2], 'W' => [0, 3],
        'K' => [2, 3], 'M' => [0, 1],
        'B' => [1, 2, 3], 'D' => [0, 2, 3],
        'H' => [0, 1, 3], 'V' => [0, 1, 2],
        'N' => ALL_STATES,
        '-' => ALL_STATES, '.' => ALL_STATES, '?' => ALL_STATES
      }.freeze
      TABLE.each_value(&:freeze)

      # Returns the Array of allowed state indices for an IUPAC
      # nucleotide symbol.  Lower case letters are accepted.
      # Raises ArgumentError for symbols that are not IUPAC
      # nucleotide codes.
      def self.state_indices(char)
        idx = TABLE[char.to_s.upcase]
        unless idx then
          raise ArgumentError,
                "unknown IUPAC nucleotide code: #{char.inspect}"
        end
        idx
      end
    end #module IUPAC

    # == Description
    #
    # Log-likelihood of a nucleotide alignment on a phylogenetic
    # tree under a substitution model (Bio::Phylo::JC69 by default),
    # with a finite mixture of rate categories.
    #
    # * The tree is a Bio::Tree with a branch length (expected number
    #   of substitutions per site) on every edge.  Multifurcations
    #   and zero-length branches are supported; a zero-length branch
    #   deterministically transmits the ancestral state (P(0) is the
    #   identity matrix), so identical bases separated only by
    #   zero-length branches behave as a single observation.
    # * The alignment maps tip labels to aligned sequence strings.
    #   Labels must match the tree leaf names exactly (one-to-one);
    #   any mismatch raises ArgumentError.
    # * IUPAC degenerate symbols (R, Y, N, ...) restrict the site to
    #   the corresponding set of concrete states.
    # * Gap policy: gaps ('-' and '.', as well as '?') are treated as
    #   fully missing data -- all four states are allowed, exactly
    #   like 'N', so the observation contributes no information.
    # * Rate categories: a finite set of relative rates, each with a
    #   weight (prior probability).  The site likelihood is the
    #   weighted average of the per-category likelihoods.  Weights
    #   are normalized to sum to 1.
    #
    # Likelihoods are computed with Felsenstein's pruning algorithm.
    # Conditional likelihood vectors are rescaled at every internal
    # node and the scaling factors are accumulated in log space, so
    # arbitrarily large trees do not underflow to -Infinity.
    #
    # For small trees, #brute_force_site_log_likelihoods explicitly
    # sums over all internal node state assignments and serves as an
    # independent cross-check of the pruning recursions.
    #
    # Because the models used here are time reversible, the position
    # of the root does not affect the likelihood (pulley principle);
    # the root is only a computational device.
    #
    # == Usage
    #
    #   lik = Bio::Phylo::Likelihood.new(tree, alignment,
    #                                    rates:   [0.5, 1.0, 2.0],
    #                                    weights: [0.25, 0.5, 0.25])
    #   lik.site_log_likelihoods   #=> Array of Float (one per column)
    #   lik.total_log_likelihood   #=> Float
    #
    class Likelihood

      # Maximum number of internal (enumerated) nodes allowed for
      # the brute force cross-check.  The number of state assignments
      # grows as 4**k, so brute force is only feasible for small
      # trees.
      BRUTE_FORCE_MAX_INTERNAL_NODES = 12

      # The Bio::Tree given to the constructor.
      attr_reader :tree

      # The substitution model (Bio::Phylo::JC69 by default).
      attr_reader :model

      # Relative rates of the rate categories (Array of Float).
      attr_reader :rates

      # Normalized weights of the rate categories (Array of Float,
      # sums to 1).
      attr_reader :weights

      # The Bio::Tree::Node used as root for the dynamic programming
      # traversal.  The likelihood itself does not depend on this
      # choice (pulley principle).
      attr_reader :root

      # Number of alignment columns (sites).
      attr_reader :number_of_sites

      # Creates a new likelihood calculator.
      #
      # tree::      a Bio::Tree whose edges all carry a branch length
      #             (Bio::Tree::Edge#distance).  Zero-length branches
      #             and multifurcations are allowed.
      # alignment:: a Hash-like object (responds to each_pair, e.g. a
      #             Hash or a Bio::Alignment with label keys) mapping
      #             tip labels (String) to aligned sequences
      #             (String).  The set of labels must be exactly the
      #             set of tree leaf names, and all sequences must
      #             have the same length.
      # model::     substitution model; must respond to
      #             state_frequencies and transition_probabilities(t).
      #             Defaults to Bio::Phylo::JC69.
      # rates::     Array of relative rates of the rate categories
      #             (each >= 0).  Defaults to [1.0] (no rate
      #             variation).
      # weights::   Array of category weights (same size as rates,
      #             each >= 0).  Normalized to sum to 1.  Defaults to
      #             uniform weights.
      # root::      Bio::Tree::Node (or a leaf/node name) to root the
      #             traversal at.  Defaults to tree.root if set,
      #             otherwise an arbitrary internal node.  The choice
      #             does not affect the likelihood.
      def initialize(tree, alignment, model: JC69.new, rates: [1.0],
                     weights: nil, root: nil)
        @tree = tree
        @model = model
        @rates = check_rates(rates)
        @weights = check_weights(weights, @rates.size)
        sequences = coerce_alignment(alignment)
        prepare_rooted_tree(root)
        assign_sequences(sequences)
        prepare_transition_matrices
      end

      # Per-site log-likelihoods: an Array of Float of length
      # number_of_sites.  A site that is impossible under the model
      # (e.g. conflicting states separated only by zero-length
      # branches) yields -Infinity.
      def site_log_likelihoods
        @site_log_likelihoods ||= compute_site_log_likelihoods
      end

      # Total log-likelihood of the alignment: the sum of the
      # per-site log-likelihoods.
      def total_log_likelihood
        site_log_likelihoods.inject(0.0) { |sum, x| sum + x }
      end
      alias log_likelihood total_log_likelihood

      # Independent cross-check of #site_log_likelihoods for small
      # trees: explicitly enumerates all assignments of states to
      # internal nodes and sums their contributions in plain (non
      # scaled) arithmetic.  Raises RuntimeError when the tree has
      # more than BRUTE_FORCE_MAX_INTERNAL_NODES internal nodes.
      def brute_force_site_log_likelihoods
        enum_nodes = @preorder.reject { |n| @leaf_states.key?(n) }
        if enum_nodes.size > BRUTE_FORCE_MAX_INTERNAL_NODES then
          raise RuntimeError,
                "brute force enumeration infeasible for " +                "#{enum_nodes.size} internal nodes " +                "(limit #{BRUTE_FORCE_MAX_INTERNAL_NODES})"
        end
        pi = @model.state_frequencies
        n_states = pi.size
        logw = @weights.map { |w| Math.log(w) }
        (0...@number_of_sites).map do |site|
          terms = @rates.each_index.map do |ci|
            lik = brute_force_site_likelihood(site, @p_matrices[ci],
                                              pi, n_states, enum_nodes)
            lik > 0.0 ? Math.log(lik) + logw[ci] : -Float::INFINITY
          end
          logsumexp(terms)
        end
      end

      # Sum of #brute_force_site_log_likelihoods.
      def brute_force_total_log_likelihood
        brute_force_site_log_likelihoods.inject(0.0) { |s, x| s + x }
      end

      private

      # Validates and freezes the rate category rates.
      def check_rates(rates)
        r = rates.to_a.map { |x| Float(x) }
        if r.empty? then
          raise ArgumentError, 'rates must not be empty'
        end
        r.each do |x|
          if x < 0.0 or x.nan? or x.infinite? then
            raise ArgumentError, "invalid rate: #{x.inspect}"
          end
        end
        r
      end

      # Validates the category weights and normalizes them to sum
      # to 1.  Without weights, uniform weights are used.
      def check_weights(weights, n_categories)
        unless weights then
          return Array.new(n_categories, 1.0 / n_categories)
        end
        w = weights.to_a.map { |x| Float(x) }
        if w.size != n_categories then
          raise ArgumentError,
                "weights size #{w.size} does not match " +                "rates size #{n_categories}"
        end
        w.each do |x|
          if x < 0.0 or x.nan? or x.infinite? then
            raise ArgumentError, "invalid weight: #{x.inspect}"
          end
        end
        sum = w.inject(0.0) { |s, x| s + x }
        unless sum > 0.0 then
          raise ArgumentError, 'weights must not all be zero'
        end
        w.map { |x| x / sum }
      end

      # Converts the alignment argument into a Hash of
      # label (String) => sequence (String) and validates it.
      def coerce_alignment(alignment)
        unless alignment.respond_to?(:each_pair) then
          raise ArgumentError,
                'alignment must map labels to sequences ' +                '(respond to each_pair)'
        end
        seqs = {}
        alignment.each_pair do |label, seq|
          label = label.to_s
          seq = seq.to_s
          if seq.empty? then
            raise ArgumentError, "empty sequence for #{label.inspect}"
          end
          if seqs.key?(label) then
            raise ArgumentError,
                  "duplicated alignment label: #{label.inspect}"
          end
          seqs[label] = seq
        end
        if seqs.empty? then
          raise ArgumentError, 'alignment is empty'
        end
        length = nil
        seqs.each do |label, seq|
          length ||= seq.length
          unless seq.length == length then
            raise ArgumentError,
                  "alignment sequences differ in length " +                  "(#{label.inspect} has #{seq.length}, " +                  "expected #{length})"
          end
        end
        seqs
      end

      # Roots the (undirected) tree and records, for every node, its
      # parent, its children and the length of the branch to its
      # parent.  Builds @preorder (parents before children) and
      # @postorder (children before parents).  The traversal is
      # iterative so that very large trees do not overflow the Ruby
      # call stack.
      def prepare_rooted_tree(root)
        nodes = @tree.nodes
        if nodes.empty? then
          raise ArgumentError, 'tree has no nodes'
        end
        @root = resolve_root(root, nodes)
        @parent = {}
        @children = Hash.new { |h, k| h[k] = [] }
        @edge_length = {}
        @parent[@root] = nil
        preorder = []
        stack = [@root]
        until stack.empty?
          node = stack.pop
          preorder << node
          @tree.adjacent_nodes(node).each do |other|
            next if other.equal?(@parent[node])
            if @parent.key?(other) then
              raise ArgumentError, 'tree has a cycle'
            end
            @parent[other] = node
            edge = @tree.get_edge(node, other)
            length = edge && edge.distance
            if length.nil? then
              raise ArgumentError,
                    "branch length missing on edge to " +                    "node #{@tree.get_node_name(other).inspect}"
            end
            length = Float(length)
            if length < 0.0 then
              raise ArgumentError,
                    "negative branch length: #{length}"
            end
            @edge_length[other] = length
            @children[node] << other
            stack.push(other)
          end
        end
        @preorder = preorder
        @postorder = preorder.reverse
      end

      # Determines the root node from the constructor argument.
      def resolve_root(root, nodes)
        case root
        when nil
          @tree.root || default_root(nodes)
        when Bio::Tree::Node
          unless @tree.include?(root) then
            raise ArgumentError, 'root node is not in the tree'
          end
          root
        else
          node = @tree.get_node_by_name(root.to_s)
          unless node then
            raise ArgumentError, "root node not found: #{root.inspect}"
          end
          node
        end
      end

      # Default root: an arbitrary internal node if any, otherwise
      # the first node of the tree.
      def default_root(nodes)
        nodes.find { |n| @tree.adjacent_nodes(n).size > 1 } || nodes.first
      end

      # Matches alignment rows to tree leaves by label and parses
      # the sequences into per-site allowed state sets.
      def assign_sequences(sequences)
        # Tips are nodes of graph degree <= 1, independently of
        # where the tree was rooted.
        tip_nodes = @preorder.select do |n|
          @tree.adjacent_nodes(n).size <= 1
        end
        by_name = {}
        tip_nodes.each do |node|
          name = @tree.get_node_name(node)
          if name.nil? or name.to_s.empty? then
            raise ArgumentError, "tree leaf has no name: #{node.inspect}"
          end
          if by_name.key?(name) then
            raise ArgumentError,
                  "duplicated leaf name in tree: #{name.inspect}"
          end
          by_name[name] = node
        end
        missing = by_name.keys - sequences.keys
        extra = sequences.keys - by_name.keys
        unless missing.empty? and extra.empty? then
          raise ArgumentError,
                "alignment/tree label mismatch: " +                "no sequence for leaves #{missing.inspect}; " +                "no leaf for sequences #{extra.inspect}"
        end
        @number_of_sites = sequences.values.first.length
        @leaf_states = {}
        by_name.each do |name, node|
          states = []
          sequences[name].each_char do |char|
            states << IUPAC.state_indices(char)
          end
          @leaf_states[node] = states
        end
      end

      # Precomputes the transition probability matrix of every
      # branch for every rate category.
      def prepare_transition_matrices
        @p_matrices = @rates.map do |rate|
          matrices = {}
          @preorder.each do |node|
            next if node.equal?(@root)
            matrices[node] =
              @model.transition_probabilities(rate * @edge_length[node])
          end
          matrices
        end
      end

      # Per-site log-likelihoods for every rate category, mixed by
      # the category weights in log space.
      def compute_site_log_likelihoods
        pi = @model.state_frequencies
        n_states = pi.size
        logs_per_category = @rates.each_index.map do |ci|
          p_matrix = @p_matrices[ci]
          (0...@number_of_sites).map do |site|
            pruning_site_log_likelihood(site, p_matrix, pi, n_states)
          end
        end
        logw = @weights.map { |w| Math.log(w) }
        (0...@number_of_sites).map do |site|
          terms = logs_per_category.each_with_index.map do |logs, ci|
            logs[site] + logw[ci]
          end
          logsumexp(terms)
        end
      end

      # Felsenstein's pruning algorithm for one site and one rate
      # category.  Conditional likelihood vectors are rescaled at
      # every internal node; the scaling factors are accumulated in
      # log space to avoid floating-point underflow on large trees.
      def pruning_site_log_likelihood(site, p_matrix, pi, n_states)
        vectors = {}
        log_scale = 0.0
        @postorder.each do |node|
          if (allowed = @leaf_states[node]) then
            vec = Array.new(n_states, 0.0)
            allowed[site].each { |s| vec[s] = 1.0 }
          else
            vec = Array.new(n_states, 1.0)
          end
          children = @children[node]
          unless children.empty? then
            children.each do |child|
              p = p_matrix[child]
              cv = vectors[child]
              n_states.times do |s|
                ps = p[s]
                sum = 0.0
                n_states.times { |sp| sum += ps[sp] * cv[sp] }
                vec[s] *= sum
              end
            end
            max = vec.max
            if max <= 0.0 then
              # Impossible data below this node (e.g. conflicting
              # states joined by zero-length branches only).
              return -Float::INFINITY
            end
            log_scale += Math.log(max)
            n_states.times { |s| vec[s] /= max }
          end
          vectors[node] = vec
        end
        root_vec = vectors[@root]
        lik = 0.0
        n_states.times { |s| lik += pi[s] * root_vec[s] }
        return -Float::INFINITY if lik <= 0.0
        Math.log(lik) + log_scale
      end

      # Likelihood of one site under one rate category, computed by
      # explicitly summing over all assignments of states to the
      # internal (non-tip) nodes.  Used as a cross-check of the
      # pruning recursions; only feasible for small trees.
      def brute_force_site_likelihood(site, p_matrix, pi, n_states,
                                      enum_nodes)
        state_of = {}
        total = 0.0
        root_states =
          if @leaf_states.key?(@root) then
            @leaf_states[@root][site]
          else
            [nil]
          end
        root_states.each do |root_state|
          state_of[@root] = root_state if root_state
          (n_states**enum_nodes.size).times do |code|
            x = code
            enum_nodes.each do |node|
              state_of[node] = x % n_states
              x /= n_states
            end
            prob = 1.0
            @preorder.each do |node|
              if node.equal?(@root) then
                prob *= pi[state_of[@root]]
              else
                ps = p_matrix[node][state_of[@parent[node]]]
                if (allowed = @leaf_states[node]) then
                  leaf_sum = 0.0
                  allowed[site].each { |s| leaf_sum += ps[s] }
                  prob *= leaf_sum
                else
                  prob *= ps[state_of[node]]
                end
              end
            end
            total += prob
          end
        end
        total
      end

      # log(sum(exp(terms))), robust to very negative terms.
      def logsumexp(terms)
        max = terms.max
        return -Float::INFINITY if max == -Float::INFINITY
        sum = 0.0
        terms.each { |t| sum += Math.exp(t - max) }
        max + Math.log(sum)
      end

    end #class Likelihood

  end #module Phylo
end #module Bio
