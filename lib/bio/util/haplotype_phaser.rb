# frozen_string_literal: true

# = bio/util/haplotype_phaser.rb - Diploid haplotype phasing by weighted MEC
#
# == Bio::HaplotypePhaser - Diploid haplotype phasing.
#
# Given sequenced fragments (reads) that each cover a set of heterozygous
# sites with allele calls 0/1 (possibly missing), find two complementary
# haplotypes and an assignment of each fragment to one of the haplotypes
# that minimizes the weighted minimum error correction (MEC) cost.

module Bio

  # = Bio::HaplotypePhaser
  #
  # Diploid haplotype phasing by the weighted minimum error correction
  # (MEC) criterion.
  #
  # Each fragment is given as a Hash mapping site indices to allele calls
  # 0/1, or as an Array of 0/1/nil calls (nil means a missing call).
  # Fragment weights must be non-negative integers (default 1).
  #
  # Sites are partitioned into connected blocks (two sites are connected
  # when some fragment covers both); disconnected blocks are solved
  # independently and the results are merged.  The exact branch-and-bound
  # search can be capped with a node budget.  When the budget is exhausted
  # before the search completes, the best solution found so far is
  # returned together with proven lower/upper bounds and Result#optimal?
  # returns false.
  #
  # == Example
  #
  #   fragments = [
  #     { 0 => 0, 1 => 0 },
  #     { 0 => 1, 1 => 1 },
  #     { 1 => 0, 2 => 1 },
  #   ]
  #   result = Bio::HaplotypePhaser.new(fragments).phase
  #   result.haplotypes   #=> e.g. [[0, 0, 1], [1, 1, nil]]
  #   result.assignments  #=> e.g. [0, 1, 0]
  #   result.cost         #=> 0
  #   result.optimal?     #=> true
  #
  class HaplotypePhaser

    # Result of a phasing run.
    class Result
      # Sorted site indices the haplotypes refer to.
      attr_reader :sites

      # Two haplotypes; arrays of 0/1/nil aligned with +sites+.
      attr_reader :haplotypes

      # Fragment -> haplotype assignment (0 or 1), aligned with the
      # input order of fragments.
      attr_reader :assignments

      # Weighted MEC cost of the returned solution.
      attr_reader :cost

      # Proven lower bound on the optimum.
      attr_reader :lower_bound

      # Upper bound on the optimum (cost of the best solution found).
      attr_reader :upper_bound

      # Per-block results (array of BlockResult), or nil.
      attr_reader :blocks

      def initialize(sites:, haplotypes:, assignments:, cost:,
                     lower_bound:, upper_bound:, complete:, blocks: nil)
        @sites = sites
        @haplotypes = haplotypes
        @assignments = assignments
        @cost = cost
        @lower_bound = lower_bound
        @upper_bound = upper_bound
        @complete = complete
        @blocks = blocks
      end

      # True when the exact search completed and the returned solution
      # is proven optimal (lower_bound == upper_bound == cost).
      def optimal?
        @complete
      end
    end

    # Result for a single connected block of sites.
    class BlockResult < Result
      # Indices of the input fragments covering the block.
      attr_reader :fragment_indices

      def initialize(fragment_indices:, **rest)
        super(**rest)
        @fragment_indices = fragment_indices
      end
    end

    # Result of an exhaustive (brute-force) verification over a small
    # fragment set.  +cost+ is the optimal weighted MEC cost and
    # +solutions+ is the number of optimal assignments, where swapping
    # the two haplotypes counts as the same solution.
    ExhaustiveResult = Struct.new(:cost, :solutions)

    # Default limit of informative fragments for .exhaustive.
    EXHAUSTIVE_MAX_FRAGMENTS = 24

    # Creates a new phaser.
    #
    # +fragments+:: Array of fragments; each fragment is a
    #               Hash[Integer => 0/1] of site calls or an Array of
    #               0/1/nil calls.
    # +weights+::   optional Array of non-negative Integer weights, one
    #               per fragment (default: all 1).
    # +budget+::    optional non-negative Integer limiting the total
    #               number of branch-and-bound search nodes across all
    #               blocks; nil (default) means unlimited.
    def initialize(fragments, weights: nil, budget: nil)
      @fragments = normalize_fragments(fragments)
      @weights = validate_weights(weights, @fragments.size)
      @budget = validate_budget(budget)
    end

    # Phases the fragments and returns a Result.
    def phase
      assignments = Array.new(@fragments.size, 0)
      budget = BudgetCounter.new(@budget)
      block_results = connected_blocks.map do |block|
        solve_block(block, budget, assignments)
      end
      merge_results(block_results, assignments)
    end

    # Exhaustively enumerates all phasings of a small fragment set and
    # returns an ExhaustiveResult with the optimal weighted MEC cost and
    # the number of optimal assignments.  Assignments that differ only by
    # swapping the two haplotypes are counted as the same solution (the
    # first informative fragment is fixed to haplotype 0).  Raises
    # ArgumentError when the number of informative fragments exceeds
    # +max_fragments+.
    def self.exhaustive(fragments, weights: nil,
                        max_fragments: EXHAUSTIVE_MAX_FRAGMENTS)
      new(fragments, weights: weights).send(:exhaustive_run, max_fragments)
    end

    private

    def exhaustive_run(max_fragments)
      informative = informative_indices
      if informative.size > max_fragments
        raise ArgumentError,
              "too many informative fragments (#{informative.size}) " \
              "for exhaustive verification (limit #{max_fragments})"
      end
      return ExhaustiveResult.new(0, 1) if informative.empty?

      best = nil
      count = 0
      (1 << (informative.size - 1)).times do |mask|
        cost = mec_cost(informative, assignment_from_mask(informative, mask))
        if best.nil? || cost < best
          best = cost
          count = 1
        elsif cost == best
          count += 1
        end
      end
      ExhaustiveResult.new(best, count)
    end

    # Assignment of the informative fragments with the first one fixed to
    # haplotype 0, decoding the remaining assignments from +mask+.
    def assignment_from_mask(informative, mask)
      assign = [0]
      (1...informative.size).each { |j| assign << ((mask >> (j - 1)) & 1) }
      assign
    end

    # Weighted MEC cost of assigning fragments +indices+ according to
    # +assign+ (0/1 per fragment).
    def mec_cost(indices, assign)
      stats = {}
      indices.each_with_index do |fi, pos|
        hap = assign[pos]
        weight = @weights[fi]
        @fragments[fi].each do |site, allele|
          col = (stats[site] ||= [[0, 0], [0, 0]])
          col[hap][allele] += weight
        end
      end
      total = 0
      stats.each_value { |col| total += col[0].min + col[1].min }
      total
    end

    def informative_indices
      (0...@fragments.size).select { |i| !@fragments[i].empty? }
    end

    def normalize_fragments(fragments)
      unless fragments.is_a?(Array)
        raise ArgumentError, 'fragments must be an Array'
      end
      fragments.map do |frag|
        case frag
        when Hash  then normalize_hash_fragment(frag)
        when Array then normalize_array_fragment(frag)
        else
          raise ArgumentError,
                "fragment must be a Hash or an Array: #{frag.inspect}"
        end
      end
    end

    def normalize_hash_fragment(frag)
      calls = {}
      frag.each do |site, allele|
        unless site.is_a?(Integer) && site >= 0
          raise ArgumentError,
                "site index must be a non-negative Integer: #{site.inspect}"
        end
        calls[site] = validate_allele(allele)
      end
      calls
    end

    def normalize_array_fragment(frag)
      calls = {}
      frag.each_with_index do |allele, site|
        calls[site] = validate_allele(allele) unless allele.nil?
      end
      calls
    end

    def validate_allele(allele)
      unless allele == 0 || allele == 1
        raise ArgumentError, "allele must be 0 or 1: #{allele.inspect}"
      end
      allele
    end

    def validate_weights(weights, size)
      return Array.new(size, 1) if weights.nil?
      unless weights.is_a?(Array) && weights.size == size
        raise ArgumentError, "weights must be an Array of size #{size}"
      end
      weights.each do |weight|
        unless weight.is_a?(Integer) && weight >= 0
          raise ArgumentError,
                "weight must be a non-negative Integer: #{weight.inspect}"
        end
      end
      weights.dup
    end

    def validate_budget(budget)
      return nil if budget.nil?
      unless budget.is_a?(Integer) && budget >= 0
        raise ArgumentError,
              "budget must be a non-negative Integer: #{budget.inspect}"
      end
      budget
    end

    # Partitions the covered sites into connected blocks.  Returns an
    # array of hashes { sites: [...], fragments: [...] } sorted by the
    # smallest site index of each block.
    def connected_blocks
      uf = UnionFind.new
      @fragments.each do |frag|
        keys = frag.keys
        keys.each { |site| uf.add(site) }
        (1...keys.size).each { |k| uf.union(keys[0], keys[k]) }
      end
      block_of = {}
      @fragments.each_with_index do |frag, i|
        next if frag.empty?
        block = (block_of[uf.find(frag.keys[0])] ||= new_block)
        block[:fragments] << i
        block[:sites].concat(frag.keys)
      end
      blocks = block_of.values
      blocks.each { |block| block[:sites] = block[:sites].uniq.sort }
      blocks.sort_by { |block| block[:sites][0] }
    end

    def new_block
      { sites: [], fragments: [] }
    end

    def solve_block(block, budget, assignments)
      sites = block[:sites]
      frag_indices = block[:fragments]
      site_pos = {}
      sites.each_with_index { |site, i| site_pos[site] = i }
      solver = BlockSolver.new(block_data(frag_indices, site_pos),
                               sites.size, budget)
      local_assign, cost, lower, complete = solver.run
      frag_indices.each_with_index do |fi, pos|
        assignments[fi] = local_assign[pos]
      end
      BlockResult.new(sites: sites, fragment_indices: frag_indices,
                      haplotypes: build_haplotypes(sites, frag_indices,
                                                   local_assign),
                      assignments: local_assign, cost: cost,
                      lower_bound: lower, upper_bound: cost,
                      complete: complete)
    end

    # Returns [frag_sites, frag_alleles, frag_weights] for the fragments
    # of a block, with site indices translated to local block positions.
    def block_data(frag_indices, site_pos)
      frag_sites = []
      frag_alleles = []
      frag_weights = []
      frag_indices.each do |fi|
        sites = []
        alleles = []
        @fragments[fi].each do |site, allele|
          sites << site_pos[site]
          alleles << allele
        end
        frag_sites << sites
        frag_alleles << alleles
        frag_weights << @weights[fi]
      end
      [frag_sites, frag_alleles, frag_weights]
    end

    # Reconstructs the two haplotypes (aligned with the block's sorted
    # sites) from a fragment assignment by weighted majority vote.
    def build_haplotypes(sites, frag_indices, local_assign)
      haps = [Array.new(sites.size), Array.new(sites.size)]
      sites.each_with_index do |site, si|
        2.times do |hap|
          w0 = 0
          w1 = 0
          frag_indices.each_with_index do |fi, pos|
            allele = local_assign[pos] == hap ? @fragments[fi][site] : nil
            next if allele.nil?
            allele.zero? ? (w0 += @weights[fi]) : (w1 += @weights[fi])
          end
          haps[hap][si] = majority_allele(w0, w1)
        end
      end
      haps
    end

    def majority_allele(w0, w1)
      if w0.zero? && w1.zero?
        nil
      elsif w1 > w0
        1
      else
        0
      end
    end

    def merge_results(block_results, assignments)
      sites = []
      hap0 = []
      hap1 = []
      cost = 0
      lower = 0
      complete = true
      block_results.each do |br|
        sites.concat(br.sites)
        hap0.concat(br.haplotypes[0])
        hap1.concat(br.haplotypes[1])
        cost += br.cost
        lower += br.lower_bound
        complete &&= br.optimal?
      end
      Result.new(sites: sites, haplotypes: [hap0, hap1],
                 assignments: assignments, cost: cost,
                 lower_bound: lower, upper_bound: cost,
                 complete: complete, blocks: block_results)
    end

    # Disjoint-set (union-find) over site indices.
    class UnionFind
      def initialize
        @parent = {}
      end

      def add(x)
        @parent[x] = x unless @parent.key?(x)
      end

      def find(x)
        root = x
        root = @parent[root] while @parent[root] != root
        until @parent[x] == root
          nxt = @parent[x]
          @parent[x] = root
          x = nxt
        end
        root
      end

      def union(a, b)
        @parent[find(a)] = find(b)
      end
    end

    # Shared counter limiting the total number of search nodes.
    class BudgetCounter
      def initialize(limit)
        @limit = limit
        @used = 0
      end

      def exhausted?
        !@limit.nil? && @used >= @limit
      end

      def consume
        @used += 1
      end
    end

    # Exact branch-and-bound MEC solver for one connected block.
    #
    # Enumerates fragment -> haplotype assignments; for a fixed
    # assignment the optimal haplotype and its cost are determined per
    # site by weighted majority, so the cost of a partial assignment is
    # a monotone lower bound used for pruning.  The first fragment is
    # fixed to haplotype 0 because swapping the two haplotypes yields an
    # equivalent solution.
    class BlockSolver
      # data:: [frag_sites, frag_alleles, frag_weights] with local site
      #        indices (see HaplotypePhaser#block_data)
      # nsites:: number of sites in the block
      # budget:: shared BudgetCounter
      def initialize(data, nsites, budget)
        frag_sites, frag_alleles, weights = data
        order = (0...frag_sites.size).sort_by { |i| [-weights[i], i] }
        @frag_sites = order.map { |i| frag_sites[i] }
        @frag_alleles = order.map { |i| frag_alleles[i] }
        @weights = order.map { |i| weights[i] }
        @order = order
        @n = order.size
        @nsites = nsites
        @budget = budget
        @w0 = Array.new(nsites) { [0, 0] }
        @w1 = Array.new(nsites) { [0, 0] }
        @assign = Array.new(@n, 0)
        @best = nil
        @best_assign = nil
        @complete = true
        @frontier_lb = nil
      end

      # Returns [assignment (original order), cost, lower_bound,
      # complete].
      def run
        search(0, 0)
        if @best_assign.nil?
          @best_assign = Array.new(@n, 0)
          @best = assignment_cost(@best_assign)
          @complete = false
          @frontier_lb = 0 if @frontier_lb.nil?
        end
        lower = @complete ? @best : [@frontier_lb, @best].min
        [unmap_assignment(@best_assign), @best, lower, @complete]
      end

      private

      def search(k, cost)
        if @budget.exhausted?
          @complete = false
          @frontier_lb = cost if @frontier_lb.nil? || cost < @frontier_lb
          return
        end
        @budget.consume
        return if @best && cost >= @best

        if k == @n
          @best = cost
          @best_assign = @assign.dup
          return
        end
        max_hap = k.zero? ? 1 : 2
        max_hap.times do |hap|
          delta = apply(k, hap, 1)
          @assign[k] = hap
          search(k + 1, cost + delta)
          apply(k, hap, -1)
        end
      end

      # Adds (sign = 1) or removes (sign = -1) fragment k on haplotype
      # +hap+.  Returns the change in total column cost when adding.
      def apply(k, hap, sign)
        weight = @weights[k] * sign
        delta = 0
        sites = @frag_sites[k]
        alleles = @frag_alleles[k]
        sites.each_with_index do |s, t|
          col0 = @w0[s]
          col1 = @w1[s]
          before = col0[hap] < col1[hap] ? col0[hap] : col1[hap]
          if alleles[t].zero?
            col0[hap] += weight
          else
            col1[hap] += weight
          end
          after = col0[hap] < col1[hap] ? col0[hap] : col1[hap]
          delta += after - before
        end
        delta * sign
      end

      # Total MEC cost of a complete assignment, computed from scratch.
      def assignment_cost(assign)
        w0 = Array.new(@nsites) { [0, 0] }
        w1 = Array.new(@nsites) { [0, 0] }
        @n.times do |k|
          hap = assign[k]
          weight = @weights[k]
          @frag_sites[k].each_with_index do |s, t|
            if @frag_alleles[k][t].zero?
              w0[s][hap] += weight
            else
              w1[s][hap] += weight
            end
          end
        end
        total = 0
        @nsites.times do |s|
          2.times do |hap|
            a = w0[s][hap]
            b = w1[s][hap]
            total += a < b ? a : b
          end
        end
        total
      end

      # Converts an assignment in solver order back to the order of the
      # fragments given to the solver.
      def unmap_assignment(assign)
        unmapped = Array.new(@n)
        @order.each_with_index { |orig, pos| unmapped[orig] = assign[pos] }
        unmapped
      end
    end
  end
end
