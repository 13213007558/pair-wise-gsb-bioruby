# frozen_string_literal: true

#
# = bio/util/rna_secondary_structure.rb - Weighted RNA secondary structure
#
# == Bio::RNASecondaryStructure
#
# Predicts an optimal pseudoknot-free RNA secondary structure under
# user-supplied constraints by maximizing a caller-configured pairing
# score. Supported constraints are forced base pairs, forbidden base
# pairs and a minimum hairpin loop length.
#
# The score is fully configured by the caller (per-pair scores). This
# class does NOT claim to implement a complete thermodynamic model.
#
# == Example
#
#  folder = Bio::RNASecondaryStructure.new("GGGAAACCC", min_loop: 3)
#  result = folder.fold
#  result.score       #=> 3 (with the default scoring of 1 per pair)
#  result.pairs       #=> [[0, 6], [1, 7], [2, 8]]
#  result.dot_bracket #=> "(((...)))"
#  result.count       #=> number of distinct optimal structures
#

module Bio
  # = Bio::RNASecondaryStructure
  #
  # Weighted Nussinov-style folding of a single RNA sequence.
  #
  # A structure is a set of base pairs (i, j) with i < j such that
  #
  # * every position is paired at most once,
  # * no two pairs cross (i.e. no pseudoknots),
  # * every hairpin loop contains at least +min_loop+ unpaired bases
  #   (that is, j - i > min_loop for every pair),
  # * every pair satisfies the caller-supplied +pair_allowed+ predicate
  #   and is not listed in +forbidden+,
  # * every pair listed in +forced+ is present.
  #
  # The score of a structure is the sum of the caller-supplied +score+
  # values of its pairs. The dynamic programming algorithm finds the
  # maximum score, one optimal structure, and the exact number of
  # distinct optimal structures. The decomposition used for counting
  # is canonical (it always branches on the leftmost position of an
  # interval), so every structure is counted exactly once, never once
  # per DP derivation.
  #
  # All positions are 0-based.
  class RNASecondaryStructure
    # Error raised when the given constraints contradict each other,
    # for example when two forced pairs cross or share a position.
    class ConstraintConflictError < ArgumentError
      # Array of human-readable descriptions of the (at least two)
      # constraints that contradict each other.
      attr_reader :constraints

      def initialize(constraints, reason = nil)
        @constraints = constraints.freeze
        message = "Conflicting constraints: #{constraints.join(' vs ')}"
        message = "#{message} (#{reason})" if reason
        super(message)
      end
    end

    # = Bio::RNASecondaryStructure::Result
    # The outcome of a folding run: one optimal structure, its score,
    # and the number of distinct optimal structures.
    class Result
      # Optimal score (Numeric), or nil when no legal structure exists.
      attr_reader :score

      # One optimal structure as a sorted array of [i, j] pairs
      # (0-based positions), or nil when infeasible.
      attr_reader :pairs

      # Dot-bracket notation of the returned structure, or nil when
      # infeasible.
      attr_reader :dot_bracket

      # Number of distinct optimal structures (Integer >= 1), or 0
      # when infeasible.
      attr_reader :count

      def initialize(score:, pairs:, count:, length:)
        @score = score
        @pairs = pairs
        @count = count
        @dot_bracket = pairs && self.class.build_dot_bracket(pairs, length)
      end

      # Returns true when at least one legal structure exists.
      def feasible?
        !@score.nil?
      end

      # Returns the dot-bracket string (or a note when infeasible).
      def to_s
        feasible? ? @dot_bracket : '(infeasible)'
      end

      # Builds the dot-bracket string for a pair list.
      def self.build_dot_bracket(pairs, length)
        chars = Array.new(length, '.')
        pairs.each do |i, j|
          chars[i] = '('
          chars[j] = ')'
        end
        chars.join
      end
    end

    # Canonical RNA base pairs (AU, GC and the GU wobble pair) used by
    # the default pairability predicate. T is treated as U.
    CANONICAL_PAIRS = {
      'A' => %w[U],
      'U' => %w[A G],
      'G' => %w[C U],
      'C' => %w[G]
    }.freeze

    # Default pairability predicate: canonical AU/GC/GU pairs.
    DEFAULT_PAIR_ALLOWED = lambda do |i, j, sequence|
      a = sequence[i].to_s.upcase
      b = sequence[j].to_s.upcase
      a = 'U' if a == 'T'
      b = 'U' if b == 'T'
      CANONICAL_PAIRS.fetch(a, []).include?(b)
    end

    # Default scoring function: 1 per pair (maximizes the number of
    # base pairs, i.e. the classic Nussinov objective).
    DEFAULT_SCORE = ->(_i, _j, _sequence) { 1 }

    # Creates a folder for the given sequence.
    #
    # === Arguments
    # * +sequence+: sequence as a String (or any object whose +to_s+
    #   gives the sequence, e.g. Bio::Sequence::NA).
    # * +min_loop+: minimum number of unpaired bases in a hairpin loop
    #   (default 3). A pair (i, j) requires j - i > min_loop.
    # * +forced+: array of [i, j] pairs (0-based) that must be present
    #   in the structure.
    # * +forbidden+: array of [i, j] pairs (0-based) that must not be
    #   present in the structure.
    # * +pair_allowed+: callable (i, j, sequence) -> truthy/falsy
    #   deciding whether positions i and j may pair. Defaults to
    #   canonical AU/GC/GU pairing.
    # * +score+: callable (i, j, sequence) -> Numeric giving the score
    #   contribution of pair (i, j). Defaults to 1 per pair.
    #
    # Raises ConstraintConflictError when the constraints contradict
    # each other, and ArgumentError for malformed constraint input.
    def initialize(sequence, min_loop: 3, forced: [], forbidden: [],
                   pair_allowed: nil, score: nil)
      @sequence = sequence.to_s
      @length = @sequence.length
      @min_loop = Integer(min_loop)
      raise ArgumentError, 'min_loop must be >= 0' if @min_loop.negative?

      @pair_allowed_fn = pair_allowed || DEFAULT_PAIR_ALLOWED
      @score_fn = score || DEFAULT_SCORE
      @allowed_cache = {}
      @score_cache = {}
      @forced = normalize_pairs(forced, 'forced')
      @forbidden = normalize_pairs(forbidden, 'forbidden')
      @forbidden_set = {}
      @forbidden.each { |a, b| @forbidden_set[[a, b]] = true }
      validate_constraints!
      @forced_partner = {}
      @forced.each do |a, b|
        @forced_partner[a] = b
        @forced_partner[b] = a
      end
      @dp_w = nil
      @dp_c = nil
    end

    # The sequence being folded (String).
    attr_reader :sequence

    # Minimum number of unpaired bases in a hairpin loop.
    attr_reader :min_loop

    # Sorted array of forced [i, j] pairs.
    attr_reader :forced

    # Sorted array of forbidden [i, j] pairs.
    attr_reader :forbidden

    # Folds the sequence and returns a Result with the optimal score,
    # one optimal structure and the number of distinct optimal
    # structures.
    def fold
      compute_dp
      if @length.zero?
        return Result.new(score: 0, pairs: [], count: 1, length: 0)
      end
      best = @dp_w[0][@length - 1]
      unless feasible_value?(best)
        return Result.new(score: nil, pairs: nil, count: 0,
                          length: @length)
      end
      pairs = []
      traceback_interval(0, @length - 1, pairs)
      Result.new(score: best, pairs: pairs.sort,
                 count: @dp_c[0][@length - 1], length: @length)
    end

    # Convenience class method: RNASecondaryStructure.fold(seq, **opts)
    # is the same as RNASecondaryStructure.new(seq, **opts).fold.
    def self.fold(sequence, **options)
      new(sequence, **options).fold
    end

    # Enumerates every legal structure (not only optimal ones) as a
    # sorted array of [i, j] pairs. Intended for short sequences: the
    # number of structures grows exponentially.
    #
    # Returns an Enumerator when no block is given.
    def each_structure
      return enum_for(:each_structure) unless block_given?

      if @length.zero?
        yield []
        return self
      end
      enumerate_interval(0, @length - 1) do |pairs|
        yield pairs.sort
      end
      self
    end

    private

    NEGATIVE_INFINITY = -Float::INFINITY

    # --- constraint validation -------------------------------------

    def normalize_pairs(pairs, kind)
      pairs.map do |pair|
        a, b = pair
        a = Integer(a)
        b = Integer(b)
        a, b = b, a if a > b
        unless a >= 0 && b < @length && a != b
          raise ArgumentError,
                "#{kind} pair (#{a}, #{b}) is out of range " \
                "for a sequence of length #{@length}"
        end
        [a, b]
      end.uniq.sort
    end

    def validate_constraints!
      @forced.each do |a, b|
        if b - a <= @min_loop
          raise ConstraintConflictError.new(
            ["forced pair (#{a}, #{b})",
             "minimum hairpin loop length #{@min_loop}"],
            'the forced pair cannot form a hairpin loop that short'
          )
        end
        unless pair_allowed?(a, b)
          raise ConstraintConflictError.new(
            ["forced pair (#{a}, #{b})", 'the base-pairing rules'],
            'the forced pair is not an allowed pair'
          )
        end
        next unless @forbidden_set[[a, b]]

        raise ConstraintConflictError.new(
          ["forced pair (#{a}, #{b})", "forbidden pair (#{a}, #{b})"],
          'the same pair is both forced and forbidden'
        )
      end
      @forced.combination(2) do |(a, b), (c, d)|
        shared = [a, b] & [c, d]
        unless shared.empty?
          raise ConstraintConflictError.new(
            ["forced pair (#{a}, #{b})", "forced pair (#{c}, #{d})"],
            "position #{shared.first} would be paired twice"
          )
        end
        next unless a < c && c < b && b < d

        raise ConstraintConflictError.new(
          ["forced pair (#{a}, #{b})", "forced pair (#{c}, #{d})"],
          'the pairs cross (would form a pseudoknot)'
        )
      end
    end

    # --- memoized caller-supplied functions -------------------------

    def pair_allowed?(i, j)
      key = i * @length + j
      if @allowed_cache.key?(key)
        @allowed_cache[key]
      else
        @allowed_cache[key] = !!@pair_allowed_fn.call(i, j, @sequence)
      end
    end

    def pair_score(i, j)
      key = i * @length + j
      if @score_cache.key?(key)
        @score_cache[key]
      else
        @score_cache[key] = @score_fn.call(i, j, @sequence)
      end
    end

    # --- dynamic programming ----------------------------------------

    # An interval [i, j] is infeasible when it contains exactly one
    # endpoint of a forced pair: the pair would have to cross the
    # interval boundary, which is impossible without a pseudoknot.
    def infeasible_interval?(i, j)
      @forced.any? do |a, b|
        (i <= a && a <= j) != (i <= b && b <= j)
      end
    end

    def feasible_value?(value)
      value != NEGATIVE_INFINITY
    end

    def w_at(i, j)
      i > j ? 0 : @dp_w[i][j]
    end

    def c_at(i, j)
      i > j ? 1 : @dp_c[i][j]
    end

    def compute_dp
      return if @dp_w

      n = @length
      @dp_w = Array.new(n) { Array.new(n, 0) }
      @dp_c = Array.new(n) { Array.new(n, 0) }
      (n - 1).downto(0) do |i|
        (i...n).each do |j|
          compute_cell(i, j)
        end
      end
    end

    # Canonical recurrence: branch on the leftmost position i of the
    # interval. Position i is either unpaired or paired with a unique
    # partner k, so every structure is generated (and counted) once.
    def compute_cell(i, j)
      if infeasible_interval?(i, j)
        @dp_w[i][j] = NEGATIVE_INFINITY
        @dp_c[i][j] = 0
        return
      end
      if (partner = @forced_partner[i])
        # i must pair with its forced partner; feasibility of the
        # interval guarantees i < partner <= j.
        @dp_w[i][j] = pair_score(i, partner) +
                      w_at(i + 1, partner - 1) + w_at(partner + 1, j)
        @dp_c[i][j] = c_at(i + 1, partner - 1) * c_at(partner + 1, j)
        return
      end
      # Option 1: i stays unpaired.
      best = w_at(i + 1, j)
      count = c_at(i + 1, j)
      # Option 2: i pairs with some k.
      (i + @min_loop + 1..j).each do |k|
        next unless pair_allowed?(i, k)
        next if @forbidden_set[[i, k]]
        next if @forced_partner[k]

        inside = w_at(i + 1, k - 1)
        outside = w_at(k + 1, j)
        next unless feasible_value?(inside) && feasible_value?(outside)

        value = pair_score(i, k) + inside + outside
        if value > best
          best = value
          count = c_at(i + 1, k - 1) * c_at(k + 1, j)
        elsif value == best
          count += c_at(i + 1, k - 1) * c_at(k + 1, j)
        end
      end
      @dp_w[i][j] = best
      @dp_c[i][j] = feasible_value?(best) ? count : 0
    end

    # --- traceback ---------------------------------------------------

    def traceback_interval(i, j, pairs)
      return if i > j

      if (partner = @forced_partner[i])
        pairs << [i, partner]
        traceback_interval(i + 1, partner - 1, pairs)
        traceback_interval(partner + 1, j, pairs)
        return
      end
      if @dp_w[i][j] == w_at(i + 1, j)
        traceback_interval(i + 1, j, pairs)
        return
      end
      (i + @min_loop + 1..j).each do |k|
        next unless pair_allowed?(i, k)
        next if @forbidden_set[[i, k]]
        next if @forced_partner[k]
        next unless pair_score(i, k) + w_at(i + 1, k - 1) +
                    w_at(k + 1, j) == @dp_w[i][j]

        pairs << [i, k]
        traceback_interval(i + 1, k - 1, pairs)
        traceback_interval(k + 1, j, pairs)
        return
      end
    end

    # --- exhaustive enumeration (short sequences) --------------------

    def enumerate_interval(i, j, &block)
      return yield [] if i > j
      return if infeasible_interval?(i, j)

      if (partner = @forced_partner[i])
        enumerate_interval(i + 1, partner - 1) do |inside|
          enumerate_interval(partner + 1, j) do |outside|
            yield [[i, partner]] + inside + outside
          end
        end
        return
      end
      # i stays unpaired.
      enumerate_interval(i + 1, j, &block)
      # i pairs with k.
      (i + @min_loop + 1..j).each do |k|
        next unless pair_allowed?(i, k)
        next if @forbidden_set[[i, k]]
        next if @forced_partner[k]

        enumerate_interval(i + 1, k - 1) do |inside|
          enumerate_interval(k + 1, j) do |outside|
            yield [[i, k]] + inside + outside
          end
        end
      end
    end
  end
end
