import React from 'react';
import { Link } from 'react-router-dom';
import Helmet from 'react-helmet';
import styles from './page2.css';
import buttonStyles from './buttons.css';


import Container from './container';

const Page2 = () => (
  <Container>
    <Helmet title="cellxgene" />
    <div className={styles.heroContainer}>
      <h1 className={styles.hero}>cellxgene</h1>
      <p className={styles.subHero}>
        A data exploration interface for single cell genetic expression matrices. Subset, filter, cluster & validate, all in one place.
      </p>
      <div>
        <Link to="/"><button className={buttonStyles.primaryButton}>explore a sample dataset</button></Link>
        <span className={styles.or}>or</span>
        <button className={buttonStyles.primaryButton}>get started with your own data</button>
      </div>
    </div>
  </Container>
);

export default Page2;
