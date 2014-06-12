package assignment4.evolution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Random;



/**
 * This is the class with the {@code main} method - it initializes the first
 * generation and then evolves it iteratively.</br>
 * To run the program, follow these steps:
 * 	Run -> Run Configurations... -> Java Application -> Right Click
 * 	-> New -> Main Class: assignment4.evolution.Evolution -> Run</br></br>
 * 
 * Your task is to implement the method {@code evolve()}.</br>
 * It has to manage the transition from one generation to the next.</br>
 * Therefore you will first need to rate the individuals (MastDistibutions)
 * based on their fitness (how well they fulfill the requirements).</br>
 * Second, you have to create a new generation from the old one by using:
 * <ul>
 * <li>selection: fit individuals survive and make it to the next generation
 * <li>crossover: features of fit individuals are recombined
 * <li>mutation: random changes occur, thus varying the genetic material
 * </ul>
 * 
 * This literature might help you get an overview of genetic algorithms or
 * evolutionary programming in general:</br>
 * 
 * <ul>
 * <li>Teske 2013: Computerprogramme züchten. In HPImgzn issue 14, pp. 50-53
 * <li>Saake, Sattler: Algorithmen und Datenstrukturen, 5th edition, pp. 86-88
 * <li>{@link http://en.wikipedia.org/wiki/Genetic_algorithm}
 * </ul>
 *
 * @author Andreas Burmeister
 * 
 * @version 0.3 05/27/14
 * 
 * @see Generation
 * @see MastDistribution
 * @see TransmitterMast
 */

public class Evolution {
	private final long budget;
	private final float minCoverage;
	private MastDistribution solution;
	private Generation currentGeneration;
	private int generationNumber;

	public Evolution(int generationSize, long budget, float minCoverage) {
		currentGeneration = new Generation(generationSize);
	    currentGeneration.randomizeAll(100); //default: 100
	    generationNumber = 1;
	    this.budget = budget;
	    this.minCoverage = minCoverage;
	    check();
	}
	
	public Generation getCurrentGeneration() {
		return currentGeneration;
	}

	public MastDistribution getSolution() {
		return solution;
	}

	public int getGenerationNumber() {
		return generationNumber;
	}
	
	private void check() {
		MastDistribution[] individuals = currentGeneration.getIndividuals();
		for (int i = 0; i < individuals.length; i++) {
			if (individuals[i].coverage() > minCoverage*CountryMap.countrySize &&
				individuals[i].getOverallCost() < budget) {
				solution = individuals[i];
				int a = solution.coverage();
				int b = CountryMap.countrySize;
				System.out.println(
					"HEUREKA!\n" +
					"Generation: " + generationNumber +
					", Individual: " + (i+1) + "\n" +
					"Coverage(3+): " + a + "/" + b + " (" +
					(int)(((float)a/b)*100) + " %)\n" +
					"Cost: " + solution.getOverallCost()
				);
				break;
			}
		}
	}
	
	public void evolve() {
		evolve(true);
	}
	
	public void evolve(boolean output) {
		///////////////////////////////////////////////////////////////////////
		//
		// TODO: implement this method (and maybe helper methods)
		//
		// How does the transition from one generation to the next happen?
		// Which individuals qualify to survive? How do you rate individuals?
		// Will mutations occur? Which kinds? Will there be crossovers? 
		//
		// e.g.:
		// crossover(currentGeneration.getIndividuals()[i], currentGeneration.getIndividuals()[j]);
		// mutate(currentGeneration.getIndividuals()[i]);
		// rateFitness(currentGeneration.getIndividuals()[i]);
		// 
		// Why do you choose to do what you do instead of something else?
		// Please comment everything you do. Help us understand your ideas.
		//
		///////////////////////////////////////////////////////////////////////
		
		// Get current and create next set of individuals
		ArrayList<MastDistribution> current = new ArrayList<MastDistribution>(Arrays.asList(currentGeneration.getIndividuals()));
		ArrayList<MastDistribution> next = new ArrayList<MastDistribution>();
		
		// Sort current generation by fitness
		Collections.sort(current, new Comparator<MastDistribution>() {
			public int compare(MastDistribution lhs, MastDistribution rhs) {
				return new Float(rateFitness(rhs)).compareTo(new Float(rateFitness(lhs)));
			}
		});
		
		// Keep some of the best without changes
		for (int i = 0; i < current.size() / 3; ++i)
			next.add(current.get(i));
		
		// Mutate all
		for (MastDistribution individual : current)
			next.add(mutate(individual));
		
		// Add some random
		for (int i = 0; i < current.size() / 2; ++i) {
			MastDistribution individual = new MastDistribution();
			individual.randomize(100);
			next.add(individual);
		}
		
		// Sort candidates for next generation by fitness
		Collections.sort(next, new Comparator<MastDistribution>() {
			public int compare(MastDistribution lhs, MastDistribution rhs) {
				return new Float(rateFitness(rhs)).compareTo(new Float(rateFitness(lhs)));
			}
		});
		
		// Print stats of best individual
		if (output) {
			MastDistribution best = next.get(0);
			float coverage = (float) best.coverage() / CountryMap.countrySize;
			float cost = (float) best.getOverallCost() / budget;
			System.out.println("Generation: " + generationNumber + ", Coverage: " + coverage +  ", Cost: " + cost);
		}
		
		// Create next generation form same amount of best individuals
		Generation nextGeneration = new Generation(current.size());
		for (int i = 0; i < current.size() && i < next.size(); ++i)
			nextGeneration.setIndividual(next.get(i), i);
		currentGeneration = nextGeneration;
		generationNumber++;
		check();
	}
	
	public float rateFitness(MastDistribution individual) {
		// Get measurements
		double cost = (double) individual.getOverallCost() / budget;
		double coverage = (double) individual.coverage() / CountryMap.countrySize;
		
		// Logarithmic cutoff for expensive distributions
		// cost = 1.0 + Math.log(Math.max(0, 1.0 - 0.9 * cost));

		// Logarithmic raise for good coverage
		// Graph at http://goo.gl/v5P9LG
		// coverage = - Math.log(Math.max(0, 1.0 - 0.9 * coverage));
		
		// Sigmoid alternative
		// Graph at http://goo.gl/n384EP
		// coverage = 1 / (1 + Math.exp(-8 * coverage + 4));
		
		// Clamp between zero and one
		// cost = clamp(cost);
		// coverage = clamp(coverage);
		
		// Weight parameters to compute fitness
		// double fitness = mix(coverage, cost, 0.05);
		// double fitness = coverage * cost;
		double fitness = coverage * Math.pow(1 - 0.15 * cost, 2.0) / (1 + cost);
		
		/*
		double fitness = 2.0 * Math.pow(coverage, 2.0)
					   + 0.9 * Math.pow(cost, 2.0)
					   + 1.3 * coverage
					   + 1.0 * cost;
		/**/
		
		/*
		// return (float) (coverage / (1 + 2.0 * cost));
		double cutoff = cost > 1.0 ? 1 / (1 + cost) : 1.0;
		coverage = 1.5 * Math.pow(coverage, 2.0) + 1.5 * coverage;
		cost = 1 + 0.5 * Math.pow(cost, 2.0) + 0.5 * cost;
		return (float) (cutoff * coverage / cost);
		/**/
		
		return (float) fitness;
	}
	
	public MastDistribution mutate(MastDistribution individual) {
		double rate = 0.1;
		
		// Get masts as hashmap
		HashMap<Point, TransmitterMast> masts = getMasts(individual);
		
		// Change amount of masts
		int amount = (int) (masts.size() * (1 + rate * random()));
		amount = Math.max(amount, 1);		
		while (amount < masts.size())
			masts.remove(pickRandomMast(masts));
		
		// Fill up with random masts
		addMasts(masts, amount - masts.size());
		
		// Randomly move masts
		final int maxMove = 10;
		int moves = (int) (masts.size() * rate * Math.random());
		int moved = 0;
		for (int abort = 0; moved < moves && abort < 10 * moves; ++abort) {
			// Grab random mast from set
			Point point = pickRandomMast(masts);
			
			// Try to find a near place
			int x = point.x + (int) random(rate * -maxMove, rate * maxMove);
			int y = point.y + (int) random(rate * -maxMove, rate * maxMove);
			Point place = new Point(x, y);
			
			if (0 > x || x > CountryMap.mapWidth - 1)
				continue;
			if (0 > y || y > CountryMap.mapHeight - 1)
				continue;
			if (CountryMap.countryMap[y][x] == 0)
				continue;
			if (masts.get(place) != null)
				continue;

			// Move mast to new place
			masts.put(place, masts.get(point));
			masts.remove(point);
			moved++;	
		}
		
		// Randomly change mast types
		// ...
		
		// Create individual from new mast set
		MastDistribution successor = new MastDistribution();
		setMasts(successor, masts);
		
		return successor;
	}
	
	public MastDistribution crossover(MastDistribution... parentIndividuals) {
		//TODO: implement
		//i.e., combine attributes from parents
		return null;
	}
	
	public boolean isTerminated() {
		//TODO: implement
		//e.g., stop after first solution (solution!=null) or try to find better ones?
		return false;
	}
	
	private HashMap<Point, TransmitterMast> getMasts(MastDistribution individual) {
		HashMap<Point, TransmitterMast> masts = new HashMap<Point, TransmitterMast>();		
		for (int y = 0; y < individual.getMastMap().length; ++y) {
			TransmitterMast[] row = individual.getMastMap()[y];
			for (int x = 0; x < row.length; ++x) {
				TransmitterMast mast = row[x];
				if (mast != null)
					masts.put(new Point(x, y), mast);
			}
		}
		return masts;
	}

	private void setMasts(MastDistribution individual, HashMap<Point, TransmitterMast> masts) {
		individual.clear();
		for (Point point : masts.keySet()) {
			try {
				individual.addMast(masts.get(point), point.x, point.y);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	private void addMasts(HashMap<Point, TransmitterMast> masts, int amount) {
	    Random random = new Random();
	    
	    // Add given amount of masts
		for (int i = 0; i < amount; ++i) {
			// Decide for mast type
			TransmitterMast mast = null;
			switch (random.nextInt(3)) {
				case 0: mast = new TransmitterMast(TransmitterMast.MastType.A); break;
				case 1: mast = new TransmitterMast(TransmitterMast.MastType.B); break;
				case 2: mast = new TransmitterMast(TransmitterMast.MastType.C); break;
			}
			
			// Find valid coordinates
			for (int abort = 0; abort < 10; ++abort) {
				int x = random.nextInt(CountryMap.mapWidth);
				int y = random.nextInt(CountryMap.mapHeight);
				Point point = new Point(x, y);
				if (CountryMap.countryMap[y][x] != 0 && masts.get(point) == null) {
					masts.put(point, mast);
					break;
				}
			};
		}
	}
	
	private Point pickRandomMast(HashMap<Point, TransmitterMast> masts) {
		Random random = new Random();
		Point[] points = masts.keySet().toArray(new Point[masts.size()]);
		return points[random.nextInt(points.length)];
	}
		
	private double random() {
		return random(-1.0f, 1.0f);
	}
	
	private double random(double min, double max) {
		return min + (max - min) * Math.random();
	}
	
	private double mix(double left, double right, double ratio) {
		return (1 - ratio) * left + ratio * right;
	}
	
	private double clamp(double value) {
		return clamp(value, 0.0, 1.0);	
	}
	
	private double clamp(double value, double min, double max) {
		return Math.min(max, Math.max(min, value));
	}

	public static void main(String[] args) {
		long budget = 100000000;
		float minCoverage = 0.9f;
		// The number of individuals in a generation. You may change this parameter
		// if you need to. However, if you choose to change it - leave a comment why.
		int generationSize = 1000;
		Evolution evolution = new Evolution(generationSize, budget, minCoverage);
		
		//stop after "terminated signal or after n iterations (default: 1000)
		while (evolution.isTerminated() && evolution.getGenerationNumber() < 1000) { 
			evolution.evolve();
    	}
	}
}