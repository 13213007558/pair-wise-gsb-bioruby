# frozen_string_literal: true

module Bio
  module DCJ
    Chromosome = Struct.new(:blocks, :circular) do
      def circular?
        circular == true
      end
    end
  end
end
