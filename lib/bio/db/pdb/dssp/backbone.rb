# frozen_string_literal: true
#
# = bio/db/pdb/dssp/backbone.rb - backbone extraction for DSSP
#
# Copyright:: Copyright (C) 2026 Naohisa Goto <ng@bioruby.org>
# License::   The Ruby License
#

require 'bio/db/pdb/dssp'

module Bio
  class PDB
    class DSSP
      # Bio::PDB::DSSP::Backbone extracts backbone coordinates from a
      # structure, detects chain breaks and reconstructs amide
      # hydrogens.  It is an internal helper of Bio::PDB::DSSP; see
      # the Bio::PDB::DSSP documentation for the exact rules.
      class Backbone
        # Backbone heavy atom names mapped to ResidueInfo accessors.
        HEAVY_ATOMS = { 'N' => :n, 'CA' => :ca, 'C' => :c, 'O' => :o }.freeze

        # Atom names recognized as backbone amide hydrogens.
        HYDROGEN_NAMES = %w[H HN].freeze

        # Assignment results grouped by chain (array of arrays of
        # ResidueInfo).
        attr_reader :chains

        # All assignment results (array of ResidueInfo).
        attr_reader :residues

        # Creates a new backbone description from _structure_ (see
        # Bio::PDB::DSSP#initialize for accepted types) using the
        # given _options_ (a hash with a :max_peptide_bond entry).
        def initialize(structure, options)
          @max_peptide_bond = options[:max_peptide_bond]
          @chains = normalize(structure).map { |residues| build_chain(residues) }
          @residues = @chains.flatten
          find_chain_breaks
          place_hydrogens
        end

        private

        def normalize(structure)
          case structure
          when Bio::PDB        then normalize_model(structure.models.first)
          when Bio::PDB::Model then normalize_model(structure)
          when Bio::PDB::Chain then [structure.residues]
          when Bio::PDB::Residue then [[structure]]
          when Array           then [structure]
          else
            raise ArgumentError, "unsupported structure: #{structure.class}"
          end
        end

        def normalize_model(model)
          return [] unless model

          model.chains.map(&:residues)
        end

        def build_chain(residues)
          residues.map { |residue| build_entry(residue) }
        end

        def build_entry(residue)
          chain = residue.chain
          entry = ResidueInfo.new(residue, chain ? chain.chain_id : nil)
          HEAVY_ATOMS.each do |name, key|
            atom = select_atom(residue, [name])
            entry.public_send(:"#{key}=", atom.xyz) if atom
          end
          entry
        end

        # Selects one atom among the given names, resolving alternate
        # locations: a blank altLoc is preferred, otherwise the
        # location with the highest occupancy wins (ties are broken
        # alphabetically).  Returns nil when no atom matches.
        def select_atom(residue, names)
          atoms = residue.atoms.select { |atom| names.include?(atom.name) }
          grouped = atoms.group_by { |atom| atom.altLoc.to_s.strip }
          best = grouped.min_by { |altloc, group| altloc_key(altloc, group) }
          best ? best.last.first : nil
        end

        def altloc_key(altloc, group)
          occupancy = group.map { |atom| atom.occupancy.to_f }.max
          [altloc.empty? ? 0 : 1, -occupancy, altloc]
        end

        def find_chain_breaks
          @chains.each do |chain|
            chain.each_cons(2) do |current, nxt|
              link(current, nxt) if peptide_bond?(current, nxt)
            end
          end
        end

        def peptide_bond?(current, nxt)
          return false unless current.c && nxt.n

          (current.c - nxt.n).r <= @max_peptide_bond
        end

        def link(current, nxt)
          current.bonded_to_next = true
          current.successor = nxt
          nxt.predecessor = current
        end

        def place_hydrogens
          @chains.each do |chain|
            chain.each_with_index do |entry, index|
              entry.h = hydrogen_position(chain, entry, index)
            end
          end
        end

        def hydrogen_position(chain, entry, index)
          return nil if entry.resname == 'PRO' || !entry.n

          atom = select_atom(entry.residue, HYDROGEN_NAMES)
          return atom.xyz if atom

          previous = index.positive? ? chain[index - 1] : nil
          return nil unless previous&.bonded_to_next && previous.c && entry.ca

          toward_c = (previous.c - entry.n).normalize
          toward_ca = (entry.ca - entry.n).normalize
          entry.n - (toward_c + toward_ca).normalize * DSSP::NH_BOND_LENGTH
        end
      end #class Backbone
    end #class DSSP
  end #class PDB
end #module Bio
