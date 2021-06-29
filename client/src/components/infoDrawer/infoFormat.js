import { Classes, H3, HTMLTable, Position, Tooltip } from "@blueprintjs/core";
import React from "react";

const ONTOLOGY_KEY = "ontology_term_id";
const COLLECTION_LINK_ORDER_BY = [
  "DOI",
  "DATA_SOURCE",
  "RAW_DATA",
  "PROTOCOL",
  "LAB_WEBSITE",
  "OTHER",
];

const buildCollectionLinks = (links) => {
  /*
  sort links by custom sort order, create view-friendly model of link types.
   */
  const sortedLinks = [...links].sort(sortCollectionLinks);
  return sortedLinks.map((link) => {
    const { link_name: name, link_type: type, link_url: url } = link;
    return {
      name: buildLinkName(name, type, url),
      type: transformLinkTypeToDisplay(type),
      url,
    };
  });
};

const buildDatasetMetadata = (singleValueCategories, corporaMetadata) => {
  /*
  transform Corpora metadata and single value categories into sort and render-friendly format.
  @returns [{key, value, tip}]
   */
  const metadata = [
    ...transformCorporaMetadata(corporaMetadata),
    ...transformSingleValueCategoriesMetadata(singleValueCategories),
  ];
  metadata.sort(sortDatasetMetadata);
  return metadata;
};

const getTableCellKeyStyles = () => {
  return { padding: "8px 8px 0 0" };
};

const getTableCellValueStyles = () => {
  return { padding: "8px 0 0 8px" };
};

const getTableStyles = () => {
  return { tableLayout: "fixed", width: "100%" };
};

const sortCollectionLinks = (l0, l1) => {
  /*
  sort collection links by custom order.
  TODO(cc) revisit - improve readability here 
   */
  return (
    COLLECTION_LINK_ORDER_BY.indexOf(l1.type) -
    COLLECTION_LINK_ORDER_BY.indexOf(l0.type)
  );
};

const buildLinkName = (name, type, url) => {
  /*
  determine name to display for collection link.
  TODO(cc) error handling
   */
  if (name) {
    return name;
  }
  if (type === "DOI") {
    return new URL(url).pathname.substring(1);
  }
  return new URL(url).host;
};

const sortDatasetMetadata = (m0, m1) => {
  /*
  sort metadata key value pairs by key - alpha, ascending
   */
  if (m0.key < m1.key) {
    return -1;
  }
  if (m0.key > m1.key) {
    return 1;
  }
  return 0;
};

const transformCorporaMetadata = (corporaMetadata) => {
  /*
  build array of view model objects from given Corpora metadata object.
  @returns [{key, value}] 
   */
  return Object.entries(corporaMetadata)
    .filter(([, value]) => {
      return value;
    })
    .map(([key, value]) => {
      return {
        key,
        value,
      };
    });
};

const transformSingleValueCategoriesMetadata = (singleValueCategories) => {
  /*
  build array of view model objects from given single value categories map, ignoring ontology terms or metadata
  without values. add ontology terms as tooltips of their corresponding values.
  @returns [{key, value, tip}] where tip is an optional ontology term for the category
   */
  return Array.from(singleValueCategories.entries())
    .filter(([key, value]) => {
      if (key.indexOf(ONTOLOGY_KEY) >= 0) {
        // skip ontology terms
        return false;
      }
      // skip metadata without values
      return value;
    })
    .map(([key, value]) => {
      const viewModel = { key, value: String(value) };
      // add ontology term as tool tip if specified
      const tip = singleValueCategories.get(`${key}_${ONTOLOGY_KEY}`);
      if (tip) {
        viewModel.tip = tip;
      }
      return viewModel;
    });
};

const transformLinkTypeToDisplay = (type) => {
  /*
  convert link type from upper snake case to title case
  TODO(cc) revisit approach here, maybe create enum-type mapping to avoid string concat inside loop?
   */
  const tokens = type.split("_");
  return tokens
    .map((token) => {
      return token.charAt(0) + token.slice(1).toLowerCase();
    })
    .join(" ");
};

const renderCollectionLinks = (collection) => {
  /*
  render collection contact and links.
  TODO(cc) handle case where there is no contact and no links?
   */
  const links = buildCollectionLinks(collection.links);
  const { contact_name: contactName, contact_email: contactEmail } = collection;
  return (
    <>
      {renderSectionTitle("Collection")}
      <HTMLTable style={getTableStyles()}>
        <tbody>
          <tr>
            <td style={getTableCellKeyStyles()}>Contact</td>
            <td style={getTableCellValueStyles()}>
              {renderCollectionContactLink(contactName, contactEmail)}
            </td>
          </tr>
          {links.map(({ name, type, url }, i) => {
            return (
              <tr {...{ key: i }}>
                <td style={getTableCellKeyStyles()}>{type}</td>
                <td style={getTableCellValueStyles()}>
                  <a href={url} rel="noopener" target="_blank">
                    {name}
                  </a>
                </td>
              </tr>
            );
          })}
        </tbody>
      </HTMLTable>
    </>
  );
};

const renderCollectionContactLink = (name, email) => {
  /*
  display collection contact's name with a link to their associated email.
   */
  if (!name && !email) {
    return null;
  }
  if (email) {
    return <a href={`mailto:${email}`}>{name}</a>;
  }
  return name;
};

const renderDatasetMetadata = (singleValueCategories, corporaMetadata) => {
  /*
  render dataset metadata, mix of meta from Corpora and attributes found in categorical field.
   */
  if (
    singleValueCategories.size === 0 &&
    Object.entries(corporaMetadata).length === 0
  ) {
    return null;
  }
  const metadata = buildDatasetMetadata(singleValueCategories, corporaMetadata);
  return (
    <>
      {renderSectionTitle("Dataset")}
      <HTMLTable style={getTableStyles()}>
        <tbody>
          {metadata.map(({ key, value, tip }) => {
            return (
              <tr {...{ key }}>
                <td style={getTableCellKeyStyles()}>{key}</td>
                <td style={getTableCellValueStyles()}>
                  <Tooltip
                    content={tip}
                    disabled={!tip}
                    minimal
                    modifiers={{ flip: { enabled: false } }}
                    position={Position.TOP}
                  >
                    {value}
                  </Tooltip>
                </td>
              </tr>
            );
          })}
        </tbody>
      </HTMLTable>
    </>
  );
};

const renderSectionTitle = (title) => {
  return (
    <p style={{ margin: "24px 0 8px" }}>
      <strong>{title}</strong>
    </p>
  );
};

const InfoFormat = React.memo(
  ({ collection, singleValueCategories, dataPortalProps = {} }) => {
    if (
      ["1.0.0", "1.1.0"].indexOf(
        dataPortalProps.version?.corpora_schema_version
      ) === -1
    ) {
      dataPortalProps = {};
    }
    const { organism } = dataPortalProps;

    return (
      <div className={Classes.DRAWER_BODY}>
        <div className={Classes.DIALOG_BODY}>
          <H3>{collection.name}</H3>
          <p>{collection.description}</p>
          {renderCollectionLinks(collection)}
          {renderDatasetMetadata(singleValueCategories, { organism })}
        </div>
      </div>
    );
  }
);

export default InfoFormat;
