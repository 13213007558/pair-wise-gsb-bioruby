# frozen_string_literal: true
#
# = bio/util/rna_fold.rb - Partition function and base-pairing probabilities
# for pseudoknot-free RNA secondary structures
#
# Copyright::   Copyright (C) 2026 BioRuby project
# License::     The Ruby License
#
# == Bio::RNAFold - Partition function and base-pairing probabilities.
#
# Instead of reporting only the minimum-energy fold, this class computes
# the partition function over all pseudoknot-free (non-crossing) secondary
# structures and derives the marginal probability of every possible base
# pair, as well as the probability of each base being unpaired.
#
# The energy model is deliberately simple: the caller provides an additive
# energy for each allowed pair (in units of kT, so that the Boltzmann
# weight of a pair is exp(-energy)) and a minimum hairpin loop length.
# The energy of a structure is the sum of its pair energies. No external
# folding program or thermodynamic parameter set is required.
#
# Forced pairs, forbidden pairs and forced-unpaired bases can be given as
# constraints. Mutually conflicting constraints are rejected by raising
# Bio::RNAFold::ConstraintError instead of silently returning all-zero
# probabilities.
#
# All internal dynamic programming is performed in log space, so the
# computed probabilities remain numerically stable even for very large
# (positive or negative) pairing energies.
#
# == Example
#
#  seq = 'GGGAAACCC'
#  # AU/GC pairs worth -2 kT, GU wobbles worth -1 kT (0-based positions)
#  energies = {}
#  seq.length.times do |i|
#    (i...seq.length).each do |j|
#      a, b = seq[i], seq[j]
#      e = case [a, b].sort.join
#          when 'AT', 'CG' then -2.0
#          when 'GT' then -1.0
#          end
#      energies[[i, j]] = e if e
#    end
#  end
#  fold = Bio::RNAFold.new(seq, pair_energy: energies, min_hairpin: 3)
#  fold.log_partition_function   #=> log of the partition function
#  fold.pair_probability(0, 8)   #=> probability that bases 0 and 8 pair
#  fold.unpaired_probability(4)  #=> probability that base 4 is unpaired
#
module Bio

  # = Bio::RNAFold
  # Partition function and base-pairing marginal probabilities for
  # pseudoknot-free RNA secondary structures with caller-supplied
  # additive pairing energies.
  class RNAFold

    # Error raised when the given constraints are mutually conflicting
    # (or admit no valid structure at all).
    class ConstraintError < ArgumentError; end

    # Default minimum number of unpaired bases inside a hairpin loop.
    DEFAULT_MIN_HAIRPIN = 3

    # Length of the sequence.
    attr_reader :length

    # Minimum number of unpaired bases inside a hairpin loop.
    attr_reader :min_hairpin

    # Forced pairs as an array of [i, j] pairs (0-based, i < j).
    attr_reader :forced_pairs

    # Forbidden pairs as an array of [i, j] pairs (0-based, i < j).
    attr_reader :forbidden_pairs

    # 0-based positions of bases that must remain unpaired.
    attr_reader :forced_unpaired

    # Creates a new Bio::RNAFold object.
    #
    # +sequence+:: the RNA sequence. Only its length and (optionally, via
    #              the energy callback) its bases are used; a String or a
    #              Bio::Sequence::NA both work.
    # +pair_energy+:: additive pairing energies in units of kT. May be
    #                 a Hash mapping [i, j] (0-based, either order) to the
    #                 energy, a 2-dimensional Array indexed as
    #                 pair_energy[i][j], or any object responding to
    #                 call(i, j). Pairs with nil, non-numeric or infinite
    #                 energy cannot be formed.
    # +min_hairpin+:: minimum number of unpaired bases in a hairpin loop;
    #                 a pair (i, j) requires j - i - 1 >= min_hairpin.
    # +forced_pairs+:: array of [i, j] pairs that must be formed.
    # +forbidden_pairs+:: array of [i, j] pairs that must not be formed.
    # +forced_unpaired+:: array of positions that must remain unpaired.
    #
    # Conflicting constraints raise Bio::RNAFold::ConstraintError.
    def initialize(sequence, pair_energy:,
                   min_hairpin: DEFAULT_MIN_HAIRPIN,
                   forced_pairs: [], forbidden_pairs: [], forced_unpaired: [])
      @sequence = sequence.to_s
      @length = @sequence.length
      @min_hairpin = Integer(min_hairpin)
      raise ArgumentError, 'min_hairpin must be non-negative' if @min_hairpin.negative?

      energy_of = normalize_energy(pair_energy)
      @forced_pairs = normalize_pairs(forced_pairs, 'forced_pairs')
      @forbidden_pairs = normalize_pairs(forbidden_pairs, 'forbidden_pairs')
      @forced_unpaired = forced_unpaired.map do |pos|
        pos = Integer(pos)
        check_position(pos, 'forced_unpaired')
        pos
      end.uniq

      build_weights(energy_of)
      validate_constraints
      apply_constraints
    end

    # Returns the caller-supplied energy (in kT units) of the pair (i, j),
    # or nil if the pair cannot be formed. Constraints are not taken into
    # account.
    def pair_energy(i, j)
      i, j = j, i if i > j
      w = @raw_log_weight[i][j]
      w.nil? || w == NEG_INFINITY ? nil : -w
    end

    # Returns the Boltzmann weight exp(-energy) of the pair (i, j) after
    # applying the constraints, or 0.0 if the pair cannot be formed.
    def pair_weight(i, j)
      i, j = j, i if i > j
      w = @log_weight[i][j]
      w.nil? || w == NEG_INFINITY ? 0.0 : Math.exp(w)
    end

    # Returns the partition function summed over all pseudoknot-free
    # structures satisfying the constraints. May be +Infinity+ for very
    # negative energies; use log_partition_function for a stable value.
    def partition_function
      Math.exp(log_partition_function)
    end

    # Returns the natural logarithm of the partition function.
    def log_partition_function
      compute
      @log_partition_function
    end

    # Returns the marginal probability that bases i and j (0-based) are
    # paired, averaged over the Boltzmann ensemble of pseudoknot-free
    # structures.
    def pair_probability(i, j)
      i, j = j, i if i > j
      check_position(i, 'i')
      check_position(j, 'j') if j != i
      compute
      @pair_prob[[i, j]] || 0.0
    end

    # Returns a Hash mapping each allowed [i, j] pair (0-based, i < j) to
    # its marginal pairing probability.
    def pair_probabilities
      compute
      @pair_prob.dup
    end

    # Returns the probability that base i (0-based) is unpaired.
    def unpaired_probability(i)
      check_position(i, 'i')
      compute
      @unpaired_prob[i]
    end

    # Returns an array with the probability of each base being unpaired.
    def unpaired_probabilities
      compute
      @unpaired_prob.dup
    end

    private

    NEG_INFINITY = -Float::INFINITY

    def normalize_energy(pair_energy)
      case pair_energy
      when Hash
        lambda do |i, j|
          pair_energy[[i, j]] || pair_energy[[j, i]]
        end
      when Array
        lambda do |i, j|
          row = pair_energy[i]
          row ? row[j] : nil
        end
      else
        unless pair_energy.respond_to?(:call)
          raise ArgumentError,
                'pair_energy must be a Hash, an Array, or respond to #call'
        end
        pair_energy
      end
    end

    def normalize_pairs(pairs, name)
      pairs.map do |pair|
        a, b = pair
        a = Integer(a)
        b = Integer(b)
        a, b = b, a if a > b
        check_position(a, name)
        check_position(b, name)
        if a == b
          raise ArgumentError, "#{name}: a base cannot pair with itself (#{a})"
        end
        [a, b]
      end.uniq
    end

    def check_position(pos, name)
      unless pos.is_a?(Integer) && pos >= 0 && pos < @length
        raise ArgumentError, "#{name}: position #{pos.inspect} out of range 0...#{@length}"
      end
    end

    def build_weights(energy_of)
      n = @length
      @raw_log_weight = Array.new(n) { Array.new(n, NEG_INFINITY) }
      n.times do |i|
        (i + @min_hairpin + 1...n).each do |j|
          energy = energy_of.call(i, j)
          next if energy.nil?
          energy = Float(energy)
          raise ArgumentError, "pair energy for (#{i}, #{j}) is NaN" if energy.nan?
          @raw_log_weight[i][j] = energy.infinite? ? NEG_INFINITY : -energy
        end
      end
      @log_weight = @raw_log_weight.map(&:dup)
    end

    def validate_constraints
      partner = {}
      @forced_pairs.each do |a, b|
        if b - a - 1 < @min_hairpin
          raise ConstraintError,
                "forced pair (#{a}, #{b}) violates min_hairpin #{@min_hairpin}"
        end
        if @raw_log_weight[a][b] == NEG_INFINITY
          raise ConstraintError, "forced pair (#{a}, #{b}) is not an allowed pair"
        end
        [a, b].each do |pos|
          if partner.key?(pos)
            raise ConstraintError,
                  "base #{pos} appears in more than one forced pair"
          end
          partner[pos] = true
        end
        if @forced_unpaired.include?(a) || @forced_unpaired.include?(b)
          raise ConstraintError,
                "forced pair (#{a}, #{b}) conflicts with forced-unpaired base"
        end
        if @forbidden_pairs.include?([a, b])
          raise ConstraintError,
                "pair (#{a}, #{b}) is both forced and forbidden"
        end
      end
      @forced_pairs.combination(2) do |pair1, pair2|
        (a, b), (c, d) = pair1[0] < pair2[0] ? [pair1, pair2] : [pair2, pair1]
        if c < b && b < d
          raise ConstraintError,
                "forced pairs (#{a}, #{b}) and (#{c}, #{d}) cross each other"
        end
      end
    end

    def apply_constraints
      n = @length
      @forbidden_pairs.each do |a, b|
        @log_weight[a][b] = NEG_INFINITY
      end
      @forced_unpaired.each do |pos|
        n.times do |k|
          @log_weight[pos][k] = NEG_INFINITY if k > pos
          @log_weight[k][pos] = NEG_INFINITY if k < pos
        end
      end
      @forced_partner = Array.new(n)
      @forced_pairs.each do |a, b|
        @forced_partner[a] = b
        @forced_partner[b] = a
        n.times do |k|
          @log_weight[a][k] = NEG_INFINITY if k > a && k != b
          @log_weight[k][a] = NEG_INFINITY if k < a
          @log_weight[b][k] = NEG_INFINITY if k > b
          @log_weight[k][b] = NEG_INFINITY if k < b && k != a
        end
      end
    end

    # Numerically stable log(exp(a) + exp(b)).
    def log_add(a, b)
      return b if a == NEG_INFINITY
      return a if b == NEG_INFINITY
      if a > b
        a + Math.log(1.0 + Math.exp(b - a))
      else
        b + Math.log(1.0 + Math.exp(a - b))
      end
    end


    # Log partition function of the complement of the interval i..j
    # (i > 0), decomposed by what happens to position i - 1.
    def outside_with_prefix(log_o, zget, i, j)
      n = @length
      partner = @forced_partner[i - 1]
      if partner
        if partner >= i && partner <= j
          # The forced partner of i - 1 lies inside the interval, so the
          # complement cannot satisfy the constraint.
          NEG_INFINITY
        elsif partner < i - 1
          @log_weight[partner][i - 1] + zget.call(partner + 1, i - 2) +
            log_o[partner][j]
        else
          @log_weight[i - 1][partner] + zget.call(j + 1, partner - 1) +
            log_o[i - 1][partner]
        end
      else
        # Position i - 1 stays unpaired.
        acc = log_o[i - 1][j]
        # Position i - 1 pairs with a position k > j, enclosing i..j.
        (j + 1...n).each do |k|
          w = @log_weight[i - 1][k]
          next if w == NEG_INFINITY
          acc = log_add(acc, w + zget.call(j + 1, k - 1) + log_o[i - 1][k])
        end
        # Position i - 1 pairs with a position k < i - 1 in the prefix.
        (0...i - 1).each do |k|
          w = @log_weight[k][i - 1]
          next if w == NEG_INFINITY
          acc = log_add(acc, w + zget.call(k + 1, i - 2) + log_o[k][j])
        end
        acc
      end
    end

    def compute
      return if @computed
      @computed = true
      n = @length
      # log_z[i][j]: log partition function of the subsequence i..j.
      log_z = Array.new(n) { Array.new(n, NEG_INFINITY) }
      zget = lambda do |i, j|
        i > j ? 0.0 : log_z[i][j]
      end
      # Inside recursion.
      (n - 1).downto(0) do |i|
        (i...n).each do |j|
          partner = @forced_partner[i]
          if partner
            log_z[i][j] =
              if partner > i && partner <= j
                @log_weight[i][partner] + zget.call(i + 1, partner - 1) +
                  zget.call(partner + 1, j)
              else
                NEG_INFINITY
              end
          else
            acc = zget.call(i + 1, j)
            (i + 1..j).each do |k|
              w = @log_weight[i][k]
              next if w == NEG_INFINITY
              acc = log_add(acc, w + zget.call(i + 1, k - 1) + zget.call(k + 1, j))
            end
            log_z[i][j] = acc
          end
        end
      end
      @log_partition_function = n.zero? ? 0.0 : log_z[0][n - 1]
      if @log_partition_function == NEG_INFINITY
        raise ConstraintError, 'constraints admit no valid structure'
      end

      # Outside recursion: log_o[i][j] is the log partition function of
      # the complement of the interval i..j, i.e. positions 0...i and
      # (j + 1)...n, where pairs may also enclose the whole interval.
      log_o = Array.new(n) { Array.new(n, NEG_INFINITY) }
      n.times do |i|
        (i...n).each do |j|
          log_o[i][j] =
            if i.zero?
              zget.call(j + 1, n - 1)
            else
              outside_with_prefix(log_o, zget, i, j)
            end
        end
      end

      # Marginal probabilities.
      log_total = @log_partition_function
      @pair_prob = {}
      paired = Array.new(n, 0.0)
      n.times do |i|
        (i + 1...n).each do |j|
          w = @log_weight[i][j]
          next if w == NEG_INFINITY
          log_p = w + zget.call(i + 1, j - 1) + log_o[i][j] - log_total
          p = log_p > 0 ? 1.0 : Math.exp(log_p)
          @pair_prob[[i, j]] = p
          paired[i] += p
          paired[j] += p
        end
      end
      @unpaired_prob = paired.map { |p| [[1.0 - p, 0.0].max, 1.0].min }
    end
  end # class RNAFold
end # module Bio
