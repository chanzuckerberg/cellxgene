import { H3, H1, UL, HTMLTable, Classes } from "@blueprintjs/core";
import React from "react";

const renderContributors = (contributors: any, affiliations: any) => {
  // eslint-disable-next-line no-constant-condition --  Temp removed contributor section to avoid publishing PII
  if (!contributors || contributors.length === 0 || true) return null;
  return (
    <>
      <H3>Contributors</H3>
      <p>
        {contributors.map((contributor: any) => {
          const { email, name, institution } = contributor;

          return (
            <span key={name}>
              {name}
              {email && `(${email})`}
              <sup>{affiliations.indexOf(institution) + 1}</sup>
            </span>
          );
        })}
      </p>
      {renderAffiliations(affiliations)}
    </>
  );
};

// generates a list of unique institutions by order of appearance in contributors
const buildAffiliations = (contributors = []) => {
  const affiliations: any = [];
  contributors.forEach((contributor) => {
    const { institution } = contributor;
    if (affiliations.indexOf(institution) === -1) {
      affiliations.push(institution);
    }
  });
  return affiliations;
};

const renderAffiliations = (affiliations: any) => {
  if (affiliations.length === 0) return null;
  return (
    <>
      <H3>Affiliations</H3>
      <UL>
        {affiliations.map((item: any, index: any) => (
          <div key={item}>
            <sup>{index + 1}</sup>
            {"  "}
            {item}
          </div>
        ))}
      </UL>
    </>
  );
};

const renderDOILink = (type: any, doi: any) => {
  if (!doi) return null;
  return (
    <>
      <H3>{type}</H3>
      <p>
        <a href={doi} target="_blank" rel="noopener">
          {doi}
        </a>
      </p>
    </>
  );
};

const ONTOLOGY_KEY = "ontology_term_id";
// Render list of metadata attributes found in categorical field
const renderDatasetMetadata = (
  singleValueCategories: any,
  corporaMetadata: any
) => {
  if (singleValueCategories.size === 0) return null;
  return (
    <>
      <H3>Dataset Metadata</H3>
      <HTMLTable
        striped
        condensed
        style={{ display: "block", width: "100%", overflowX: "auto" }}
      >
        <thead>
          <tr>
            <th>Field</th>
            <th>Label</th>
            <th>Ontology ID</th>
          </tr>
        </thead>
        <tbody>
          {Object.entries(corporaMetadata).map(([key, value]) => {
            return (
              <tr {...{ key }}>
                <td>{`${key}:`}</td>
                {/* @ts-expect-error ts-migrate(2322) FIXME: Type 'unknown' is not assignable to type 'ReactNod... Remove this comment to see the full error message */}
                <td>{value}</td>
                <td />
              </tr>
            );
          })}
          {Array.from(singleValueCategories).reduce((elems, pair) => {
            // @ts-expect-error ts-migrate(2488) FIXME: Type 'unknown' must have a '[Symbol.iterator]()' m... Remove this comment to see the full error message
            const [category, value] = pair;
            // If the value is empty skip it
            if (!value) return elems;

            // If this category is a ontology term, let's add its value to the previous node
            if (String(category).includes(ONTOLOGY_KEY)) {
              const prevElem = (elems as any).pop();
              const newChildren = [...prevElem.props.children];
              newChildren.splice(2, 1, [<td key="ontology">{value}</td>]);
              // Props aren't extensible so we must clone and alter the component to append the new child
              (elems as any).push(
                React.cloneElement(prevElem, prevElem.props, newChildren)
              );
            } else {
              // Create the list item
              (elems as any).push(
                <tr key={category}>
                  <td>{`${category}:`}</td>
                  <td>{value}</td>
                  <td />
                </tr>
              );
            }
            return elems;
          }, [])}
        </tbody>
      </HTMLTable>
    </>
  );
};

// Renders any links found in the config where link_type is not "SUMMARY"
// If there are no links in the config, render the aboutURL
const renderLinks = (projectLinks: any, aboutURL: any) => {
  if (!projectLinks && !aboutURL) return null;
  if (projectLinks)
    return (
      <>
        <H3>Project Links</H3>
        <UL>
          {projectLinks.map((link: any) => {
            if (link.link_type === "SUMMARY") return null;
            return (
              <li key={link.link_name}>
                <a href={link.link_url} target="_blank" rel="noopener">
                  {link.link_name}
                </a>
              </li>
            );
          })}
        </UL>
      </>
    );

  return (
    <>
      <H3>More Info</H3>
      <p>
        <a href={aboutURL} target="_blank" rel="noopener">
          {aboutURL}
        </a>
      </p>
    </>
  );
};

const InfoFormat = React.memo(
  // @ts-expect-error ts-migrate(2339) FIXME: Property 'datasetTitle' does not exist on type '{ ... Remove this comment to see the full error message
  ({ datasetTitle, singleValueCategories, aboutURL, dataPortalProps = {} }) => {
    if (
      ["1.0.0", "1.1.0"].indexOf(
        dataPortalProps.version?.corpora_schema_version
      ) === -1
    ) {
      dataPortalProps = {};
    }
    const {
      title,
      publication_doi: doi,
      preprint_doi: preprintDOI,
      organism,
      contributors,
      project_links: projectLinks,
    } = dataPortalProps;

    const affiliations = buildAffiliations(contributors);

    return (
      <div className={Classes.DIALOG_BODY}>
        <div className={Classes.DIALOG_BODY}>
          <H1>{title ?? datasetTitle}</H1>
          {renderContributors(contributors, affiliations)}
          {renderDatasetMetadata(singleValueCategories, { organism })}
          {renderLinks(projectLinks, aboutURL)}
          {renderDOILink("DOI", doi)}
          {renderDOILink("Preprint DOI", preprintDOI)}
        </div>
      </div>
    );
  }
);

export default InfoFormat;
