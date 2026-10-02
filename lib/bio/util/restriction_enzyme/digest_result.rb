#
# bio/util/restriction_enzyme/digest_result.rb - Result of explicit digestion
#

# Copyright:: Copyright (C) 2026 BioRuby project
# License:: Ruby License
#

require 'bio/util/restriction_enzyme'

module Bio

class RestrictionEnzyme

# Sticky or blunt end produced by one paired double-stranded cut.
StickyEnd = Struct.new(:kind, :length, :sequence) do
  def blunt?
    kind == :blunt
  end

  def sticky?
    !blunt?
  end
end

# A paired cut in sequence boundary coordinates.
DigestCut = Struct.new(:enzyme, :orientation, :primary_cut, :complement_cut,
                       :recognition_start, :recognition_end, :sticky_end)

# A recognition site present in the sequence but blocked by a caller rule.
BlockedSite = Struct.new(:enzyme, :orientation, :start, :end,
                         :methylations)

# A double-stranded fragment. Strand strings are padded with spaces so their
# base pairs align in the same column.
DigestFragment = Struct.new(:primary, :complement, :topology,
                            :primary_start, :primary_end,
                            :complement_start, :complement_end,
                            :left_cut, :right_cut) do
  def circular?
    topology == :circular
  end

  def linear?
    topology == :linear
  end

  def left_sticky_end
    left_cut ? left_cut.sticky_end : nil
  end

  def right_sticky_end
    right_cut ? right_cut.sticky_end : nil
  end
end

# Result returned by RestrictionEnzyme.digest.
class DigestResult
  attr_reader :sequence, :circular, :cuts, :fragments, :sticky_ends,
              :blocked_sites

  def initialize(sequence, circular)
    @sequence = sequence
    @circular = circular
    @cuts = []
    @fragments = []
    @sticky_ends = []
    @blocked_sites = []
  end

  def circular?
    @circular
  end

  alias blocked_recognition_sites blocked_sites
end

end # RestrictionEnzyme
end # Bio
