#
# bio/util/restriction_enzyme/digest_rule.rb - Explicit restriction enzyme rule
#
# Copyright:: Copyright (C) 2026 BioRuby project
# License:: Ruby License
#

require 'bio/sequence'
require 'bio/util/restriction_enzyme'

module Bio

class RestrictionEnzyme

# A caller supplied restriction enzyme rule used by RestrictionEnzyme.digest.
#
# Cuts are represented in 0-based boundary coordinates: a cut at index n is
# between bases n and n+1. Negative indices and indices at or beyond the
# recognition site are used for cuts outside the recognition site.
class DigestRule
  # Human-readable name used in result records.
  attr_reader :name

  # Recognition site in 5' to 3' orientation on the primary strand.
  attr_reader :recognition

  # Boundary cuts local to the aligned rule: [[primary, complement], ...].
  attr_reader :cuts

  # Methylation positions that block cleavage.
  # Entries are [strand, position, modification].
  attr_reader :sensitive_to

  # +strand+ is +:both+ (default), +:+,
  # or :- and controls whether the reverse complement is scanned.
  attr_reader :strand

  def self.cast(rule)
    return rule if rule.is_a?(DigestRule)

    case rule
    when Hash
      pattern = rule.fetch(:recognition) { rule.fetch(:pattern) }
      if pattern.to_s.include?('^')
        cast_rule = from_pattern(pattern, rule[:name])
        cast_rule.configure(rule[:cuts] || rule[:cut],
                            rule[:sensitive_to] || rule[:sensitivity],
                            rule[:strand])
        cast_rule
      else
        new(pattern,
            name: rule[:name],
            cuts: rule[:cuts] || rule[:cut],
            sensitive_to: rule[:sensitive_to] || rule[:sensitivity],
            strand: rule[:strand])
      end
    when DoubleStranded
      from_double_stranded(rule)
    else
      from_pattern(rule)
    end
  end

  def self.from_pattern(pattern, name = nil)
    enzyme = DoubleStranded.new(pattern)
    rule = from_double_stranded(enzyme)
    rule.instance_variable_set(:@name, name || pattern)
    rule
  end

  def self.from_double_stranded(enzyme)
    aligned = enzyme.aligned_strands
    primary_start = aligned.primary.index(enzyme.primary.stripped)
    recognition = enzyme.primary.stripped.to_s
    cuts = enzyme.cut_locations.map do |cut|
      [cut.primary - primary_start,
       cut.complement - primary_start]
    end
    new(recognition, name: nil, cuts: cuts)
  end

  def initialize(recognition,
                 name: nil, cuts: nil, sensitive_to: nil, strand: :both)
    @recognition = Bio::Sequence::NA.new(recognition).to_s
    raise ArgumentError, 'recognition site is empty' if @recognition.empty?

    @name = name
    @cuts = normalize_cuts(cuts)
    @sensitive_to = normalize_sensitivity(sensitive_to)
    @strand = (strand || :both).to_sym
    unless [:both, :+, :-].include?(@strand)
      raise ArgumentError, "unknown strand: #{@strand.inspect}"
    end
  end

  def scan_plus?
    @strand == :both || @strand == :+
  end

  def scan_minus?
    @strand == :both || @strand == :-
  end

  def palindromic?
    reverse = Bio::Sequence::NA.new(recognition).reverse_complement.to_s
    @palindromic ||= recognition == reverse
  end

  def label
    @name || recognition
  end

  def configure(cuts = nil, sensitivity = nil, strand = nil)
    @cuts = normalize_cuts(cuts) unless cuts.nil?
    @sensitive_to = normalize_sensitivity(sensitivity) unless sensitivity.nil?
    @strand = (strand || @strand).to_sym
    self
  end

  def blocking_methylations(methylations, orientation, position, size,
                            circular = false)
    sensitive_to.flat_map do |strand, local_position, modification|
      physical_strand = physical_strand(orientation, strand)
      physical_position = position + local_strand_position(orientation,
                                                           local_position,
                                                           size)
      physical_position %= size if circular
      methylations.select do |methylation|
        methylation.strand == physical_strand &&
          methylation.position == physical_position &&
          modification_matches?(modification, methylation.modification)
      end
    end
  end

  def blocked_by?(methylations, orientation, position, size)
    !blocking_methylations(methylations, orientation, position, size).empty?
  end

  def physical_cuts(position, orientation)
    cuts.map do |primary_cut, complement_cut|
      if orientation == :+
        [position + primary_cut, position + complement_cut]
      else
        [position + minimum_boundary + maximum_boundary - primary_cut,
         position + minimum_boundary + maximum_boundary - complement_cut]
      end
    end
  end

  def maximum_boundary
    @maximum_boundary ||= [recognition.length - 1, *cuts.flatten].max
  end

  def minimum_boundary
    @minimum_boundary ||= [0, *cuts.flatten].min
  end

  protected

  def physical_strand(orientation, strand)
    if orientation == :+ || strand.nil?
      strand
    else
      strand == :+ ? :- :(:+)
    end
  end

  def modification_matches?(expected, actual)
    expected.nil? || expected.to_s.casecmp(actual.to_s).zero?
  end

  def local_strand_position(orientation, local_position, size)
    if orientation == :+
      local_position
    else
      size - 1 - local_position
    end
  end

  def normalize_cuts(cuts)
    cuts = default_cuts if cuts.nil?
    cuts = [cuts] unless cuts.first.is_a?(Array)
    cuts.map do |cut|
      cut = Array(cut)
      unless cut.size == 2
        raise ArgumentError, 'a cut requires two strand boundaries'
      end
      cut.map(&:to_i)
    end
  end

  def default_cuts
    [[recognition.length / 2 - 1, recognition.length / 2 - 1]]
  end

  def normalize_sensitivity(sensitivity)
    Array(sensitivity).map do |entry|
      case entry
      when Hash
        [entry.fetch(:strand).to_sym,
         entry.fetch(:position).to_i,
         entry[:modification] || entry[:type]]
      when Array
        [entry[0].to_sym, entry[1].to_i, entry[2]]
      else
        raise ArgumentError, "invalid sensitivity rule: #{entry.inspect}"
      end
    end
  end
end

# Methylation on a specific strand and 0-based base position.
Methylation = Struct.new(:strand, :position, :modification)

end # RestrictionEnzyme
end # Bio
