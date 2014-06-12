package assignment4.evolution;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

/**
 * A class that opens a frame in order to show the {@code generationSize}
 * individuals of the {@code currentGeneration} and has a button to trigger the
 * {@code evolve} method.</br>
 * There are a couple of parameters for the visuals you may change.</br>
 * 
 * To run the program, follow these steps:
 * 	Run -> Run Configurations... -> Java Application -> Right Click
 * 	-> New -> Main Class: assignment4.evolution.EvolutionGUI -> Run</br></br>
 *
 * @author Andreas Burmeister
 * 
 * @version 0.6 05/27/14
 * 
 * @see Evolution
 */

public class EvolutionGUI implements ActionListener {
	final Evolution evolution;
	private JFrame f = new JFrame("Evolution");
	private JButton nextButton;
	

	public EvolutionGUI(Evolution evolution) {
		this.evolution = evolution;
	}
	
	public void showGUI() {
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		EvoPanel evoPanel = new EvoPanel();
		nextButton = new JButton("next generation");
		nextButton.setActionCommand("evolve");
		nextButton.addActionListener(this);
		evoPanel.add(nextButton);
		f.add(evoPanel);
		f.pack();
		f.setResizable(true);
	    f.setVisible(true);
	}
	
	public void actionPerformed(ActionEvent e) {
		if("evolve".equals(e.getActionCommand()) && !evolution.isTerminated()) {
			evolution.evolve();
			// disable button if terminated
			nextButton.setEnabled(!evolution.isTerminated());
			f.repaint();
		}
		
		if("evolve".equals(e.getActionCommand()) && !evolution.isTerminated()) {
			
		}
	}


	@SuppressWarnings("serial")
	private class EvoPanel extends JPanel {
		
		public EvoPanel() {
			setBorder(BorderFactory.createLineBorder(Color.black));
		}
		
		
		/* You may change the parameters BELOW. **********************************/
		
		/**
		 * @param size: size of a single area
		 * @param padding: distance between two areas
		 * @param spacing: border around areas on map
		 * @param border: border around maps in frame
		 */
		private int size = 6;     //default:  7
		private int padding = 2;  //default:  2
		private int spacing = 10; //default: 10
		private int border = 5;   //default:  5
		
		/**
		 * @param xNumber: number of maps in a row
		 */
		private int xNumber = 6;  //default:  4
		
		/**
		 * @param colors: colors for the areas inside the country
		 * @param color: color for areas outside the country
		 */
		private int[][] colors = {
		    // r,   g,   b
			{255,   0,   0}, // coverage 0 , default: 255,   0,   0 (red)
			{255, 127,   0}, // coverage 1 , default: 255, 127,   0 (orange)
			{255, 255,   0}, // coverage 2 , default: 255, 255,   0 (yellow)
			{ 63, 255,   0}, // coverage 3 , default:  63, 255,   0 (light green)
			{  0, 191,   0}  // coverage 4+, default:   0, 191,   0 (green)
		};
		private int[] color = {
			// r,   g,   b
			 191, 191, 191 //default: 191, 191, 191
		};
		
		/* You may change the parameters ABOVE. **********************************/
		
		
		private int yNumber = (int)(Math.ceil(evolution.getCurrentGeneration().size()*1D/xNumber));
		//number of columns needed
		
		private int xSize = //horizontal size of a single map
			CountryMap.mapWidth*size+(CountryMap.mapWidth-1)*padding+2*spacing;
		private int ySize = //vertical size of a single map
			CountryMap.mapHeight*size+(CountryMap.mapHeight-1)*padding+2*spacing;
		
		public Dimension getPreferredSize() {
			return new Dimension(
				xNumber*xSize+(xNumber+1)*border,
				yNumber*ySize+(yNumber+1)*border);
		}
		
		public void paintComponent(Graphics g) {
		//Draws the maps for the current generation of individuals.
			super.paintComponent(g);
			for (int i = 0; i < evolution.getCurrentGeneration().size(); i++) {
				MastDistribution individual = evolution.getCurrentGeneration().getIndividuals()[i];
				int xCoord = //
					((i%xNumber)+1)*border+
					(i%xNumber)*xSize;
				int yCoord = //
					((int)Math.ceil(i/xNumber)+1)*border+
					(int)Math.floor(i/xNumber)*ySize;
				g.setColor(Color.white);
				g.fillRect(xCoord, yCoord, xSize, ySize);
				g.setColor(Color.black);
				g.drawRect(xCoord, yCoord, xSize, ySize);
				int[][] coveragemap = individual.getCoverageMap();
				int[][] countrymap =
					CountryMap.countryMap;
				for(int row = 0; row < CountryMap.mapHeight; row++) {
					for(int col = 0; col < CountryMap.mapWidth; col++) {
						int coverage =
							Math.min(coveragemap[row][col], this.colors.length-1);
						int xPos = col*(size+padding)+spacing+xCoord;
						int yPos = row*(size+padding)+spacing+yCoord;
						if (countrymap[row][col] == 0) {
							g.setColor(new Color(
								this.color[0],
								this.color[1],
								this.color[2]));
						} else {
							g.setColor(new Color(
								this.colors[coverage][0],
								this.colors[coverage][1],
								this.colors[coverage][2]));
						}
						g.fillRect(xPos, yPos, size, size);
						if (countrymap[row][col] == 1) {
							g.setColor(Color.black);
							g.drawRect(xPos, yPos, size, size);
						}
					}
				}
				g.setColor(Color.white);
				g.fillRect(xCoord+spacing, yCoord+spacing, 6*size+5*padding, 2*size+padding);
				g.setColor(Color.black);
				g.drawRect(xCoord+spacing, yCoord+spacing, 6*size+5*padding, 2*size+padding);
				g.drawString(String.format("%.5f", evolution.rateFitness(individual)), 
						xCoord+padding+spacing, yCoord+padding+2*spacing);
			}
		}
	}
	
	public static void main(String[] args) {
		long budget = 100000000;
		float minCoverage = 0.9f;
		// The number of individuals in a generation. You may change this parameter
		// if you need to - but be warned that the GUI will not be able to display
		// many individuals and will only become slower.
		// However, if you choose to change it - leave a comment why.
		//default: 18
		int generationSize = 18;
		Evolution evolution = new Evolution(generationSize, budget, minCoverage);
		final EvolutionGUI gui = new EvolutionGUI(evolution);
		
		// Give some evolutions in advance
		int runs = 10000;
 		for (int i = 0; i < runs - 1; ++i)
 			evolution.evolve(false);
 		System.out.println("Finished " + runs + " runs in advance.");
		evolution.evolve(true);
 		
		// Click your way through the generations until the evolution is terminated
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				gui.showGUI();
			}
		});
	}
}