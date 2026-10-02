# frozen_string_literal: true

#
# = bio/util/rna_secondary_structure.rb - Pseudoknot-free RNA secondary
#   structure prediction under a weighted base-pairing model
#
# Copyright::   Copyright (C) 2026 BioRuby Project
# License::     The Ruby License
#

module Bio
  # = Bio::RNASecondaryStructure
  #
  # Predicts an optimal pseudoknot-free RNA secondary structure under a
  # weighted base-pairing model, given experimentally derived constraints.
  #
  # The score of a structure is simply the sum of the weights of its base
  # pairs.  The weights are fully configured by the caller (see the block
  # form of ::new); this class does NOT implement a thermodynamic energy
  # model.  Unpaired bases contribute nothing to the score.
  #
  # Supported constraints:
  #
  # * <tt>min_hairpin:</tt> minimum number of unpaired bases in a hairpin
  #   loop.  A pair (i, j) is allowed only when <tt>j - i - 1 >=
  #   min_hairpin</tt>.
  # * <tt>forced:</tt> pairs that must appear in the structure (e.g.
  #   experimentally confirmed pairings).
  # * <tt>forbidden:</tt> pairs that must not appear in the structure.
  #
  # Forced pairs that contradict each other (crossing pairs, pairs sharing
  # a base) or contradict another constraint (a forbidden pair, the
  # minimum hairpin length) are rejected up front with an ArgumentError
  # that names the two conflicting constraints.
  #
  # Besides one optimal structure, the result reports how many *distinct*
  # optimal structures exist.  The dynamic programming recursion decides
  # the fate of the leftmost position of each interval (unpaired, or
  # paired with exactly one partner), so every structure corresponds to
  # exactly one decomposition and is counted exactly once.
  #
  # == Example
  #
  #   folder = Bio::RNASecondaryStructure.new("GGGAAACCC", min_hairpin: 3)
  #   result = folder.fold
  #   result.score       #=> 3
  #   result.count       #=> 1
  #   result.pairs       #=> [[0, 8], [1, 7], [2, 6]]
  #   result.dot_bracket #=> "(((...)))"
  #
  #   # caller-configured scoring
  #   weights = { "GC" => 3, "CG" => 3, "AU" => 2, "UA" => 2,
  #               "GU" => 1, "UG" => 1 }
  #   seq = "GGGAAACCC"
  #   result = Bio::RNASecondaryStructure.fold(seq) do |i, j|
  #     weights[seq[i, 1] + seq[j, 1]]
  #   end
  #
  class RNASecondaryStructure
    # Negative infinity used for infeasible (sub)problems.
    NEGATIVE_INFINITY = -Float::INFINITY

    # Default pair weights: canonical Watson-Crick pairs and G-U wobble
    # pairs score 1, everything else is not allowed to pair.
    DEFAULT_PAIR_WEIGHTS = {
      'AU' => 1, 'UA' => 1,
      'GC' => 1, 'CG' => 1,
      'GU' => 1, 'UG' => 1
    }.freeze

    # The outcome of a folding computation.
    #
    # When the constraints are infeasible (no structure satisfies them),
    # #feasible? is false, #count is 0, #score is -Float::INFINITY and
    # both #pairs and #dot_bracket are nil.
    class Result
      # Score of the optimal structure(s).
      attr_reader :score
      # Number of distinct optimal structures.
      attr_reader :count
      # One optimal structure, as a sorted array of [i, j] base pairs
      # (zero-based positions, i < j).  nil when infeasible.
      attr_reader :pairs

      def initialize(score, count, pairs, length) # :nodoc:
        @score = score
        @count = count
        @pairs = pairs
        @length = length
      end

      # Returns true when at least one structure satisfies the constraints.
      def feasible?
        @count.positive?
      end

      # Dot-bracket notation of the optimal structure in #pairs, e.g.
      # "(((...)))".  Returns nil when infeasible.
      def dot_bracket
        return nil unless feasible?

        chars = Array.new(@length, '.')
        @pairs.each do |i, j|
          chars[i] = '('
          chars[j] = ')'
        end
        chars.join
      end
    end

    # The RNA sequence as an upcased String (T is treated as U).
    attr_reader :sequence
    # Minimum number of unpaired bases in a hairpin loop.
    attr_reader :min_hairpin
    # Forced pairs, as sorted [[i, j], ...] with zero-based positions.
    attr_reader :forced_pairs
    # Forbidden pairs, as sorted [[i, j], ...] with zero-based positions.
    attr_reader :forbidden_pairs

    # Creates a folder for the given sequence.
    #
    # _sequence_:: the RNA sequence (String or any object responding to
    #              to_s, e.g. Bio::Sequence::NA).
    # <tt>min_hairpin:</tt>:: minimum number of unpaired bases in a
    #                         hairpin loop (default 3).
    # <tt>forced:</tt>:: array of [i, j] pairs (zero-based) that must
    #                    appear in the structure.
    # <tt>forbidden:</tt>:: array of [i, j] pairs (zero-based) that must
    #                       not appear in the structure.
    #
    # The optional block configures the scoring: it is called with two
    # zero-based positions (i, j) and must return the numeric weight of
    # the pair, or nil (or false) when the pair is not allowed.  When no
    # block is given, DEFAULT_PAIR_WEIGHTS is used.
    #
    # Raises ArgumentError when constraints contradict each other; the
    # message names the two conflicting constraints.
    def initialize(sequence, min_hairpin: 3, forced: [], forbidden: [], &weight)
      @sequence = sequence.to_s.upcase.tr('T', 'U')
      @length = @sequence.length
      @min_hairpin = validate_min_hairpin(min_hairpin)
      @weight = weight || default_weight_proc
      @forbidden_pairs = normalize_pairs(forbidden, 'forbidden')
      @forbidden_set = @forbidden_pairs.each_with_object({}) do |pair, hash|
        hash[pair] = true
      end
      @forced_pairs = normalize_pairs(forced, 'forced')
      validate_forced_pairs!
      @forced_partner = {}
      @forced_pairs.each do |a, b|
        @forced_partner[a] = b
        @forced_partner[b] = a
      end
    end

    # Folds the sequence and returns a Result.
    def fold
      @result ||= compute_result
    end

    # Convenience shortcut: Bio::RNASecondaryStructure.fold(seq, **opts) { }
    # is equivalent to Bio::RNASecondaryStructure.new(seq, **opts) { }.fold.
    def self.fold(sequence, **options, &weight)
      new(sequence, **options, &weight).fold
    end

    private

    def validate_min_hairpin(min_hairpin)
      value = Integer(min_hairpin)
      raise ArgumentError, "min_hairpin must be >= 0: #{value}" if value.negative?

      value
    rescue TypeError
      raise ArgumentError, "min_hairpin must be an Integer: #{min_hairpin.inspect}"
    end

    def normalize_pairs(pairs, kind)
      pairs.map { |pair| normalize_pair(pair, kind) }.uniq
    end

    def normalize_pair(pair, kind)
      a, b = convert_pair(pair, kind)
      a, b = b, a if a > b
      if a == b || a.negative? || b >= @length
        raise ArgumentError,
              "invalid #{kind} pair (#{a}, #{b}) for a sequence of length #{@length}"
      end
      [a, b]
    end

    def convert_pair(pair, kind)
      a, b = pair
      [Integer(a), Integer(b)]
    rescue TypeError, ArgumentError
      raise ArgumentError, "invalid #{kind} pair: #{pair.inspect}"
    end

    def validate_forced_pairs!
      @forced_pairs.each { |pair| validate_forced_pair!(*pair) }
      validate_forced_pairs_compatible!
    end

    def validate_forced_pair!(a, b)
      loop_size = b - a - 1
      if loop_size < @min_hairpin
        raise ArgumentError,
              "forced pair (#{a}, #{b}) conflicts with minimum hairpin length " +              "#{@min_hairpin}: the hairpin loop would contain only #{loop_size} base(s)"
      end
      return unless @forbidden_set[[a, b]]

      raise ArgumentError,
            "forced pair (#{a}, #{b}) conflicts with forbidden pair (#{a}, #{b})"
    end

    def validate_forced_pairs_compatible!
      used = {}
      @forced_pairs.each do |a, b|
        [a, b].each do |pos|
          other = used[pos]
          if other
            raise ArgumentError,
                  "forced pair (#{a}, #{b}) conflicts with forced pair " +                  "(#{other[0]}, #{other[1]}): position #{pos} is used by both"
          end
          used[pos] = [a, b]
        end
      end
      @forced_pairs.sort.combination(2) do |(a, b), (c, d)|
        next unless a < c && c < b && b < d

        raise ArgumentError,
              "forced pair (#{a}, #{b}) conflicts with forced pair (#{c}, #{d}): " +              'the pairs cross, which would form a pseudoknot'
      end
    end

    def default_weight_proc
      sequence = @sequence
      lambda do |i, j|
        DEFAULT_PAIR_WEIGHTS[sequence[i, 1] + sequence[j, 1]]
      end
    end

    def compute_result
      return Result.new(0, 1, [], 0) if @length.zero?

      compute_tables
      score = @best[0][@length - 1]
      ways = @ways[0][@length - 1]
      pairs = ways.positive? ? traceback_pairs : nil
      Result.new(score, ways, pairs, @length)
    end

    def compute_tables
      @best = Array.new(@length) { Array.new(@length, NEGATIVE_INFINITY) }
      @ways = Array.new(@length) { Array.new(@length, 0) }
      (@length - 1).downto(0) do |i|
        (i...@length).each { |j| compute_cell(i, j) }
      end
    end

    # Fills @best[i][j] and @ways[i][j]: the best score achievable on the
    # interval i..j and the number of distinct structures achieving it.
    def compute_cell(i, j)
      partner = @forced_partner[i]
      if partner
        compute_forced_cell(i, j, partner)
        return
      end
      # Option 1: position i stays unpaired.
      best_score, best_ways = value_of(i + 1, j)
      # Option 2: position i pairs with some k.
      ((i + @min_hairpin + 1)..j).each do |k|
        option = pair_option(i, j, k)
        next unless option

        score, ways = option
        if score > best_score
          best_score = score
          best_ways = ways
        elsif score == best_score
          best_ways += ways
        end
      end
      @best[i][j] = best_score
      @ways[i][j] = best_ways
    end

    # A forced position i must pair with its forced partner; if the
    # partner lies outside the interval, the subproblem is infeasible.
    def compute_forced_cell(i, j, partner)
      option = partner < i || partner > j ? nil : pair_option(i, j, partner)
      if option
        @best[i][j], @ways[i][j] = option
      else
        @best[i][j] = NEGATIVE_INFINITY
        @ways[i][j] = 0
      end
    end

    # Score and count for interval i..j; the empty interval scores 0 and
    # has exactly one (empty) structure.
    def value_of(i, j)
      return [0, 1] if i > j

      [@best[i][j], @ways[i][j]]
    end

    # Returns [score, ways] for placing pair (i, k) inside interval i..j,
    # or nil when the pair is not allowed here.
    def pair_option(i, j, k)
      return nil if @forbidden_set[[i, k]]

      other = @forced_partner[k]
      return nil if other && other != i

      weight = @weight.call(i, k)
      return nil if weight.nil? || weight == false
      return nil if weight.is_a?(Float) && !weight.finite?

      left_score, left_ways = value_of(i + 1, k - 1)
      return nil if left_ways.zero?

      right_score, right_ways = value_of(k + 1, j)
      return nil if right_ways.zero?

      [weight + left_score + right_score, left_ways * right_ways]
    end

    def traceback_pairs
      pairs = []
      trace(0, @length - 1, pairs)
      pairs.sort
    end

    # Reconstructs one optimal structure on interval i..j.  Only called
    # on feasible intervals.
    def trace(i, j, pairs)
      return if i > j

      partner = @forced_partner[i]
      if partner
        pairs << [i, partner]
        trace(i + 1, partner - 1, pairs)
        trace(partner + 1, j, pairs)
        return
      end
      target = @best[i][j]
      left_score, = value_of(i + 1, j)
      if left_score == target
        trace(i + 1, j, pairs)
        return
      end
      trace_pair_options(i, j, target, pairs)
    end

    def trace_pair_options(i, j, target, pairs)
      ((i + @min_hairpin + 1)..j).each do |k|
        option = pair_option(i, j, k)
        next unless option

        score, = option
        next unless score == target

        pairs << [i, k]
        trace(i + 1, k - 1, pairs)
        trace(k + 1, j, pairs)
        return
      end
    end
  end
end
