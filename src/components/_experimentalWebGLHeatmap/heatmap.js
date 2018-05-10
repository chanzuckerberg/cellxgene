// jshint esversion: 6
import React from "react";
import _ from "lodash";
import DeckGL, {
  PointCloudLayer,
  ScreenGridLayer,
  COORDINATE_SYSTEM
} from "deck.gl";
import OrbitController from "./orbit-control";
import { Popup } from "./popup";

class Heatmap extends React.Component {
  constructor(props) {
    super(props);

    this.onChangeViewport = this.onChangeViewport.bind(this);
    this.onInitialized = this.onInitialized.bind(this);
    this.onResize = this.onResize.bind(this);
    this.onUpdate = this.onUpdate.bind(this);
    this.onHover = this.onHover.bind(this);

    this.state = {
      width: 0,
      height: 0,
      points: [],
      sampleExpressionMatrix: [
        { color: [0, 255, 0], position: [100, 100] },
        { color: [0, 255, 0], position: [100, 100] },
        { color: [0, 255, 0], position: [100, 100] }
      ],
      progress: 0,
      popup: {
        displayed: false,
        x: 0,
        y: 0,
        title: ""
      },
      viewport: {
        lookAt: [0, 0, 0],
        distance: 1,
        rotationX: 0,
        rotationY: 0,
        fov: 30,
        minDistance: 0.5,
        maxDistance: 3
      }
    };
  }

  getColor(cluster) {
    let color = [0, 0, 0];

    switch (cluster) {
      case 0:
        color = [166, 206, 227];
        break;
      case 1:
        color = [31, 120, 180];
        break;
      case 2:
        color = [178, 223, 138];
        break;
      case 3:
        color = [51, 160, 44];
        break;
      case 4:
        color = [251, 154, 153];
        break;
      case 5:
        color = [227, 26, 28];
        break;
      case 6:
        color = [253, 191, 111];
        break;
      case 7:
        color = [255, 127, 0];
        break;
      case 8:
        color = [202, 178, 214];
        break;
      case 9:
        color = [106, 61, 154];
        break;
    }

    return color;
  }

  componentWillMount() {
    window.addEventListener("resize", this.onResize);
    this.onResize();
  }

  componentDidMount() {
    this.canvas.fitBounds([-0.5, -0.5, -0.5], [0.5, 0.5, 0.5]);

    this.fetchData();

    window.requestAnimationFrame(this.onUpdate);
  }

  componentWillUnmount() {
    window.removeEventListener("resize", this.onResize);
  }

  fetchData() {
    fetch(
      "https://raw.githubusercontent.com/zdenekhynek/data-science-capstone-visualisation/master/public/clusters.json"
    ).then(res => {
      res.json().then(obj => {
        const clusters = Object.keys(obj).map(k => obj[k]);
        const points = clusters.map(cluster => {
          const position = [cluster.x, cluster.y, cluster.z];
          const color = [255, 0, 0];
          const id = cluster.id;
          const title = cluster.webTitle;
          return { id, title, position, color };
        });

        this.setState({ points, progress: 1 });
      });
    });
  }

  onHover(d) {
    let popup = { displayed: false };
    console.log("d", d);

    if (d.object) {
      const object = d.object;
      popup = {
        id: object.id,
        title: object.title,
        displayed: true,
        x: d.x,
        y: d.y
      };
    }

    this.setState({ popup });
  }

  onResize() {
    const { innerWidth: width, innerHeight: height } = window;
    this.setState({ width: width / 1.5, height: height / 1.5 });
  }

  onInitialized(gl) {
    gl.clearColor(0, 0, 0, 1);
    gl.enable(gl.DEPTH_TEST);
    gl.depthFunc(gl.LEQUAL);
  }

  onChangeViewport(viewport) {
    this.setState({
      rotating: !viewport.isDragging,
      viewport: { ...this.state.viewport, ...viewport }
    });
  }

  onUpdate() {
    const { viewport } = this.state;
    window.requestAnimationFrame(this.onUpdate);
  }

  renderPointCloudLayer() {
    return (
      this.state.points.length &&
      new PointCloudLayer({
        id: "point-cloud-layer",
        data: this.state.points,
        projectionMode: COORDINATE_SYSTEM.IDENTITY,
        pickable: true,
        onHover: this.onHover,
        getPosition: d => d.position,
        getNormal: d => [0, 0.5, 0.2],
        getColor: d => d.color,
        radiusPixels: 2
      })
    );
  }

  renderGridLayer() {
    /**
     * Data format:
     * [
     *   {position: [-122.4, 37.7]},
     *   ...
     * ]
     */
    const screenGridLayer = new ScreenGridLayer({
      id: "screen-grid-layer",
      data: this.state.sampleExpressionMatrix,
      projectionMode: COORDINATE_SYSTEM.IDENTITY,
      pickable: true,
      getPosition: d => d.position,
      getColor: d => d.color,
      cellSizePixels: 40
    });

    return screenGridLayer;
  }

  renderDeckGLCanvas() {
    const { width, height, viewport } = this.state;
    const canvasProps = { width, height, ...viewport };
    const glViewport = OrbitController.getViewport(canvasProps);

    return (
      width &&
      height && (
        <OrbitController
          {...canvasProps}
          ref={canvas => {
            this.canvas = canvas;
          }}
          onChangeViewport={this.onChangeViewport}
        >
          <DeckGL
            width={width}
            height={height}
            viewport={glViewport}
            layers={[
              // this.renderPointCloudLayer(),
              this.renderGridLayer()
            ].filter(Boolean)}
            onWebGLInitialized={this.onInitialized}
          />
        </OrbitController>
      )
    );
  }

  render() {
    const { width, height, popup } = this.state;
    if (!width || !height) {
      return null;
    }

    const renderedPopup = popup.displayed ? <Popup {...popup} /> : null;

    return (
      <div id="heatmap">
        {this.renderDeckGLCanvas()}
        {renderedPopup}
      </div>
    );
  }
}

export default Heatmap;
