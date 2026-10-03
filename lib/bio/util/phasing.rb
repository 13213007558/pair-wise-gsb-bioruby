# frozen_string_literal: true
#
# = bio/util/phasing.rb - Diploid haplotype phasing by minimum error correction
#
# Copyright::   Copyright (C) 2026 BioRuby project
# License::     The Ruby License
#
# == Bio::Phasing - diploid haplotype phasing by weighted minimum error
# correction (MEC).
#
# A set of sequencing fragments each covers a number of heterozygous
# sites with allele calls 0/1 (a call may also be missing). The goal is
# to find two complementary haplotypes and an assignment of every
# fragment to one of the haplotypes that minimizes the weighted error
# correction cost: the sum, over all fragments, of the fragment weight
# times the number of called sites at which the fragment disagrees with
# its assigned haplotype.
#
# The problem is solved exactly by branch-and-bound over fragment
# assignments. Per-column majority voting is not used, because it would
# lose the association of alleles along each fragment. Sites that are
# not connected by any fragment form independent blocks and are solved
# separately. An exact-search budget (maximum number of search nodes)
# can be given; when the budget is exhausted before optimality is
# proven, the best solution found so far is returned together with a
# proven lower bound, and Result#optimal? returns false.
#
# == Example
#
#  fragments = [
#    { 0 => 0, 1 => 0, 2 => 1 },
#    { 1 => 0, 2 => 1 },
#    { 0 => 1, 1 => 1, 2 => 0 },
#  ]
#  result = Bio::Phasing.new(fragments).phase
#  result.optimal?     #=> true
#  result.cost         #=> 0
#  result.haplotypes   #=> [[0, 0, 1], [1, 1, 0]]
#  result.assignment   #=> [0, 0, 1]
#

module Bio

  # = Bio::Phasing
  #
  # Diploid haplotype phasing by weighted minimum error correction
  # (MEC). See the documentation of bio/util/phasing.rb for an
  # overview.
  #
  class Phasing

    # Default exact-search budget: the maximum total number of
    # branch-and-bound search nodes explored by #phase.
    DEFAULT_MAX_NODES = 1_000_000

    # Result of a phasing run returned by Bio::Phasing#phase.
    class Result
      # Array of BlockResult, one per connected block of sites.
      attr_reader :blocks
      # Best (smallest) weighted error correction cost found. This is
      # an upper bound of the optimal cost.
      attr_reader :cost
      # Proven lower bound of the optimal cost.
      attr_reader :lower_bound
      # Total number of branch-and-bound search nodes explored.
      attr_reader :nodes

      def initialize(blocks, n_fragments) # :nodoc:
        @blocks = blocks
        @n_fragments = n_fragments
        @cost = blocks.inject(0) { |sum, b| sum + b.cost }
        @lower_bound = blocks.inject(0) { |sum, b| sum + b.lower_bound }
        @nodes = blocks.inject(0) { |sum, b| sum + b.nodes }
      end

      # Returns true when the returned solution is proven to be
      # optimal, i.e. the exact search completed within the node
      # budget. Returns false when the budget was exhausted; in that
      # case #cost is the best solution found so far (an upper bound)
      # and #lower_bound is a proven lower bound of the optimum.
      def optimal?
        @blocks.all?(&:optimal?)
      end

      # Returns the sorted Array of all covered heterozygous site
      # indices, across all blocks.
      def sites
        @sites ||= @blocks.flat_map(&:sites).sort
      end

      # Returns the two phased haplotypes as [h1, h2], each an Array
      # of 0/1 alleles aligned with #sites. h2 is the complement of
      # h1. The relative orientation of disconnected blocks is
      # arbitrary.
      def haplotypes
        @haplotypes ||= build_haplotypes
      end

      # Returns an Array of 0/1 haplotype assignments, one per input
      # fragment (0 = haplotype 1, 1 = haplotype 2). Fragments without
      # any called site are assigned 0.
      def assignment
        return @assignment if @assignment
        @assignment = Array.new(@n_fragments, 0)
        @blocks.each do |block|
          block.assignment.each { |i, side| @assignment[i] = side }
        end
        @assignment
      end

      private

      def build_haplotypes
        alleles = [{}, {}]
        @blocks.each do |block|
          h1, h2 = block.haplotypes
          block.sites.each_with_index do |site, i|
            alleles[0][site] = h1[i]
            alleles[1][site] = h2[i]
          end
        end
        [0, 1].map { |k| sites.map { |site| alleles[k][site] } }
      end
    end

    # Phasing result of one connected block of heterozygous sites.
    class BlockResult
      # Sorted Array of site indices covered by this block.
      attr_reader :sites
      # Two haplotypes [h1, h2] as Arrays of 0/1 aligned with #sites.
      attr_reader :haplotypes
      # Hash mapping input fragment index to a 0/1 assignment.
      attr_reader :assignment
      # Best weighted error correction cost found for this block.
      attr_reader :cost
      # Proven lower bound of the optimal cost of this block.
      attr_reader :lower_bound
      # Number of branch-and-bound nodes explored for this block.
      attr_reader :nodes

      def initialize(attributes) # :nodoc:
        @sites = attributes.fetch(:sites)
        @haplotypes = attributes.fetch(:haplotypes)
        @assignment = attributes.fetch(:assignment)
        @cost = attributes.fetch(:cost)
        @lower_bound = attributes.fetch(:lower_bound)
        @nodes = attributes.fetch(:nodes)
        @optimal = attributes.fetch(:optimal)
      end

      # Returns true when this block's solution is proven optimal.
      def optimal?
        @optimal
      end
    end

    # Result of Bio::Phasing#exhaustive.
    class ExhaustiveResult
      # Optimal weighted error correction cost.
      attr_reader :cost
      # Number of optimal solutions. Solutions that differ only by
      # swapping the two haplotypes (within each block) are the same
      # solution and are counted only once.
      attr_reader :optimal_count

      def initialize(cost, optimal_count) # :nodoc:
        @cost = cost
        @optimal_count = optimal_count
      end
    end

    # Normalized input fragment. For internal use only.
    class Fragment # :nodoc:
      attr_reader :index, :weight, :sites, :alleles

      def initialize(index, weight, sites, alleles)
        @index = index
        @weight = weight
        @sites = sites
        @alleles = alleles
      end
    end

    # Creates a phasing problem.
    #
    # +fragments+:: Array of fragment call maps. Each call map is a
    #               Hash (or an Array of [site, allele] pairs) mapping
    #               an Integer site index to an allele call 0 or 1.
    #               Sites absent from the map (or mapped to nil) are
    #               treated as missing (uncalled).
    # +weights+::   optional Array of non-negative Integer fragment
    #               weights. Every fragment has weight 1 when omitted.
    def initialize(fragments, weights: nil)
      validate_input_sizes!(fragments, weights)
      @fragments = fragments.each_with_index.map do |calls, i|
        weight = weights ? weights[i] : 1
        validate_weight!(weight)
        sites, alleles = normalize_calls(calls)
        Fragment.new(i, weight, sites, alleles)
      end
    end

    # Solves the phasing problem and returns a Result.
    #
    # +max_nodes+ is the exact-search budget: the maximum total number
    # of branch-and-bound nodes explored over all blocks
    # (DEFAULT_MAX_NODES by default). When nil or Float::INFINITY, the
    # search is unbounded and the solution is always proven optimal.
    # When the budget is exhausted first, the best solution found so
    # far is returned and Result#optimal? is false.
    def phase(max_nodes: DEFAULT_MAX_NODES)
      remaining = check_budget(max_nodes)
      blocks = fragment_blocks.map do |group|
        solver = BlockSolver.new(group, remaining)
        result = solver.solve
        remaining = [remaining - solver.nodes, 0].max
        result
      end
      Result.new(blocks, @fragments.length)
    end

    # Solves the problem by exhaustive enumeration of all fragment
    # assignments and returns an ExhaustiveResult. Intended for
    # verifying optimal values on small fragment sets. Raises
    # ArgumentError when a connected block has more than
    # +max_fragments+ fragments.
    def exhaustive(max_fragments: 20)
      cost = 0
      count = 1
      fragment_blocks.each do |group|
        check_block_size!(group, max_fragments)
        block_cost, block_count = exhaustive_block(group)
        cost += block_cost
        count *= block_count
      end
      ExhaustiveResult.new(cost, count)
    end

    # Same as Bio::Phasing.new(fragments, weights: weights).exhaustive.
    def self.exhaustive(fragments, weights: nil, max_fragments: 20)
      new(fragments, weights: weights).exhaustive(max_fragments: max_fragments)
    end

    private

    def validate_input_sizes!(fragments, weights)
      fragments.is_a?(Array) or
        raise ArgumentError, 'fragments must be an Array'
      return if weights.nil?
      weights.is_a?(Array) && weights.length == fragments.length or
        raise ArgumentError,
              'weights must be an Array of the same length as fragments'
    end

    def validate_weight!(weight)
      return if weight.is_a?(Integer) && weight >= 0
      raise ArgumentError,
            "fragment weight must be a non-negative Integer: #{weight.inspect}"
    end

    def check_budget(max_nodes)
      return Float::INFINITY if max_nodes.nil?
      unless max_nodes.is_a?(Numeric) && max_nodes >= 0
        raise ArgumentError,
              "max_nodes must be a non-negative number or nil: #{max_nodes.inspect}"
      end
      max_nodes
    end

    def check_block_size!(group, max_fragments)
      return if group.length <= max_fragments
      raise ArgumentError,
            "block of #{group.length} fragments exceeds " \
            "max_fragments=#{max_fragments}"
    end

    def normalize_calls(calls)
      pairs = {}
      each_call(calls) { |site, allele| add_call(pairs, site, allele) }
      sorted = pairs.sort
      [sorted.map(&:first), sorted.map(&:last)]
    end

    def each_call(calls)
      return calls.each_pair { |s, a| yield s, a } if calls.is_a?(Hash)
      unless calls.is_a?(Array)
        raise ArgumentError,
              "fragment calls must be a Hash or an Array: #{calls.inspect}"
      end
      calls.each { |pair| yield valid_pair(pair) }
    end

    def valid_pair(pair)
      pair.is_a?(Array) && pair.length == 2 or
        raise ArgumentError, "invalid [site, allele] pair: #{pair.inspect}"
      pair
    end

    def add_call(pairs, site, allele)
      site.is_a?(Integer) or
        raise ArgumentError, "site index must be an Integer: #{site.inspect}"
      return if allele.nil? # missing call
      allele == 0 || allele == 1 or
        raise ArgumentError, "allele call must be 0 or 1: #{allele.inspect}"
      pairs[site] = allele
    end

    # Groups fragments into connected blocks: two fragments are in the
    # same block when they are linked through shared sites. Fragments
    # without any called site belong to no block.
    def fragment_blocks
      parent = build_site_unions
      groups = Hash.new { |hash, key| hash[key] = [] }
      @fragments.each do |f|
        groups[uf_find(parent, f.sites.first)] << f unless f.sites.empty?
      end
      groups.values
    end

    def build_site_unions
      parent = {}
      @fragments.each do |f|
        next if f.sites.empty?
        root = uf_find(parent, f.sites.first)
        f.sites.each { |site| uf_union(parent, root, site) }
      end
      parent
    end

    def uf_find(parent, x)
      root = x
      root = parent[root] while parent.key?(root)
      y = x
      while parent.key?(y)
        parent[y], y = root, parent[y]
      end
      root
    end

    def uf_union(parent, root, site)
      site_root = uf_find(parent, site)
      parent[site_root] = root unless site_root == root
    end

    # Enumerates all 2**(n-1) assignments of a block (the first
    # fragment is fixed to side 0, because swapping the two haplotypes
    # gives the same solution) and returns [cost, number_of_optima].
    def exhaustive_block(group)
      frags = group.sort_by(&:index)
      best = nil
      count = 0
      (1 << (frags.length - 1)).times do |mask|
        cost = block_cost(frags, mask_sides(frags.length, mask))
        best, count = update_optimum(best, count, cost)
      end
      [best, count]
    end

    def mask_sides(n, mask)
      sides = [0]
      (n - 1).times { |j| sides << ((mask >> j) & 1) }
      sides
    end

    def update_optimum(best, count, cost)
      return [cost, 1] if best.nil? || cost < best
      return [best, count + 1] if cost == best
      [best, count]
    end

    def block_cost(frags, sides)
      counts = Hash.new { |hash, key| hash[key] = [[0, 0], [0, 0]] }
      frags.each_with_index do |f, i|
        f.sites.each_with_index do |site, k|
          counts[site][sides[i]][f.alleles[k]] += f.weight
        end
      end
      counts.values.inject(0) { |sum, (c0, c1)| sum + c0.min + c1.min }
    end

    # Exact branch-and-bound solver for one connected block of
    # fragments. For internal use only.
    class BlockSolver # :nodoc:
      # Number of search nodes explored.
      attr_reader :nodes

      def initialize(fragments, max_nodes)
        @frags = fragments.sort_by { |f| [f.sites.first, f.sites.last, f.index] }
        @n = @frags.length
        @max_nodes = max_nodes
        index_sites
        reset_search_state
      end

      # Runs the search and returns a BlockResult.
      def solve
        dfs(0)
        assign = @best_assign || Array.new(@n, 0)
        cost = @best_cost || assignment_cost(assign)
        lower = @cut ? [@frontier_min || 0, cost].min : cost
        BlockResult.new(sites: @sites, haplotypes: build_haplotypes(assign),
                        assignment: index_assignment(assign), cost: cost,
                        lower_bound: lower, nodes: @nodes, optimal: !@cut)
      end

      private

      def index_sites
        @site_frags = Hash.new { |hash, key| hash[key] = [] }
        @frags.each_with_index do |f, i|
          f.sites.each_with_index do |site, k|
            @site_frags[site] << [i, f.alleles[k]]
          end
        end
        @sites = @site_frags.keys.sort
        @counts = {}
        @sites.each { |site| @counts[site] = [[0, 0], [0, 0]] }
      end

      def reset_search_state
        @lb = 0
        @assign = Array.new(@n)
        @best_cost = nil
        @best_assign = nil
        @nodes = 0
        @cut = false
        @frontier_min = nil
      end

      # Depth-first search over fragment assignments. The lower bound
      # of a partial assignment is the sum, over all sites and both
      # sides, of the smaller of the two allele weights; it can only
      # grow as more fragments are assigned.
      def dfs(i)
        return update_best if i == @n
        child_sides(i).each do |side|
          place(i, side)
          visit_child(i)
          unplace(i)
        end
      end

      def update_best
        return if @best_cost && @lb >= @best_cost
        @best_cost = @lb
        @best_assign = @assign.dup
      end

      def child_sides(i)
        return [0] if i.zero? # haplotype-swap symmetry
        [0, 1].sort_by { |side| [delta_for(i, side), side] }
      end

      def visit_child(i)
        child_lb = @lb
        return note_frontier(child_lb) if @cut || child_lb >= best_so_far
        return cut_search(child_lb) if @nodes >= @max_nodes
        @nodes += 1
        dfs(i + 1)
      end

      def cut_search(child_lb)
        @cut = true
        note_frontier(child_lb)
      end

      def best_so_far
        @best_cost || Float::INFINITY
      end

      def note_frontier(lb)
        @frontier_min = lb if @frontier_min.nil? || lb < @frontier_min
      end

      def with_counts(i, side)
        f = @frags[i]
        f.sites.each_with_index do |site, k|
          yield @counts[site][side], f.alleles[k], f.weight
        end
      end

      def place(i, side)
        delta = 0
        with_counts(i, side) do |c, allele, w|
          before = c.min
          c[allele] += w
          delta += c.min - before
        end
        @lb += delta
        @assign[i] = side
      end

      def unplace(i)
        delta = 0
        with_counts(i, @assign[i]) do |c, allele, w|
          before = c.min
          c[allele] -= w
          delta += c.min - before
        end
        @lb += delta
        @assign[i] = nil
      end

      def delta_for(i, side)
        delta = 0
        with_counts(i, side) do |c, allele, w|
          before = c.min
          c[allele] += w
          delta += c.min - before
          c[allele] -= w
        end
        delta
      end

      def assignment_cost(assign)
        @sites.inject(0) { |sum, site| sum + site_cost(site, assign) }
      end

      def site_cost(site, assign)
        c = [[0, 0], [0, 0]]
        @site_frags[site].each do |pos, allele|
          c[assign[pos]][allele] += @frags[pos].weight
        end
        c[0].min + c[1].min
      end

      def build_haplotypes(assign)
        h1 = @sites.map { |site| majority_allele(site, assign) }
        [h1, h1.map { |a| 1 - a }]
      end

      def majority_allele(site, assign)
        w0 = 0
        w1 = 0
        @site_frags[site].each do |pos, allele|
          next unless assign[pos].zero?
          if allele.zero?
            w0 += @frags[pos].weight
          else
            w1 += @frags[pos].weight
          end
        end
        w0 >= w1 ? 0 : 1
      end

      def index_assignment(assign)
        assignment = {}
        @frags.each_with_index { |f, i| assignment[f.index] = assign[i] }
        assignment
      end
    end
  end
end
