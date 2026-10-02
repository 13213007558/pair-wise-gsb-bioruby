# frozen_string_literal: true

# = bio/dcj.rb - Double-cut-and-join distance and scenarios

require 'set'
require_relative 'dcj/endpoint'
require_relative 'dcj/chromosome'
require_relative 'dcj/operation'
require_relative 'dcj/genome_parser'
require_relative 'dcj/adjacencies'
require_relative 'dcj/chromosome_traversal'
require_relative 'dcj/chromosome_normalizer'
require_relative 'dcj/chromosome_builder'
require_relative 'dcj/operation_validator'
require_relative 'dcj/operation_candidates'
require_relative 'dcj/distance_calculator'
require_relative 'dcj/scenario_builder'
require_relative 'dcj/genome'

module Bio
  module DCJ
    def self.distance(source, target)
      source = Genome.new(source) unless source.is_a?(Genome)
      target = Genome.new(target) unless target.is_a?(Genome)
      validate_pair(source, target)
      DistanceCalculator.new(source, target).call
    end

    def self.scenario(source, target)
      source = Genome.new(source) unless source.is_a?(Genome)
      target = Genome.new(target) unless target.is_a?(Genome)
      validate_pair(source, target)
      ScenarioBuilder.new(source, target).call
    end

    def self.validate_pair(source, target)
      return if source.block_ids == target.block_ids

      raise ArgumentError, 'genomes must contain the same gene blocks'
    end
  end
end
